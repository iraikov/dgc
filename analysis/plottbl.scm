(require-extension matchable)
(require-library srfi-1 srfi-4 irregex data-structures files posix extras ploticus)
(import
 (only srfi-4 list-f64vector)
 (only srfi-1 filter list-tabulate)
 (only files make-pathname)
 (only posix glob)
 (only data-structures ->string alist-ref compose)
 (only extras fprintf random)
 (only mathh cosh tanh log10)
 (prefix ploticus plot:)
 )


(define comment-pat (string->irregex "^%.*"))

(define (sample n v)
  (let ((ub (vector-length v)))
    (let ((idxs (list-tabulate n (lambda (i) (random ub)))))
      (map (lambda (i) (vector-ref v i)) idxs)
      ))
    )


(define (fresults1 port)
  (lambda (vs)
    (for-each
     (lambda (v)
       (let (
             (input-resistance (f64vector-ref v 1))
             (membrane-tau (f64vector-ref v 6))
             (spike-threshold (f64vector-ref v 9))
             (spike-amplitude(f64vector-ref v 10))
             (spike-ahp (f64vector-ref v 11))
             (rel-amplitude-dend1 (f64vector-ref v 12))
             (rel-amplitude-dend2 (f64vector-ref v 13))
             (rel-amplitude-dend3 (f64vector-ref v 14))
             (rel-amplitude-dend4 (f64vector-ref v 15))
             (rel-amplitude-dend5 (f64vector-ref v 16))
             (number-of-spikes (f64vector-ref v 17))
             (mean-firing-rate (f64vector-ref v 18))
             (mean-isi (f64vector-ref v 19))
             (isi-adaptation3 (f64vector-ref v 23))
             (isi-adaptation4 (f64vector-ref v 24))
             )
         (fprintf port
                  "~A,~A,~A,~A,~A,~A,~A,~A,~A,~A,~A,~A,~A,~A,~A~%" 
                  input-resistance
                  membrane-tau
                  spike-threshold
                  spike-amplitude
                  spike-ahp
                  rel-amplitude-dend1
                  rel-amplitude-dend2
                  rel-amplitude-dend3
                  rel-amplitude-dend4
                  rel-amplitude-dend5
                  number-of-spikes
                  mean-firing-rate
                  mean-isi
                  isi-adaptation3
                  isi-adaptation4)
         
         ))
     vs)
    ))

(define (tbl-plot plot-label . data-files)

  (let 

   (
    (data
     (fold
      (lambda (dat-file data)
        (let ((data1 (map (lambda (line) 
                            (list->f64vector (map string->number (string-split line ","))))
                          (filter (lambda (line) (not (irregex-match comment-pat line)))
                                  (read-lines dat-file)))))
          (cons data1 data)))
      '()
      data-files))
    )

  (let-values (
               ((fd1 temp-path1) (file-mkstemp "/tmp/tbl-plot.s1.XXXXXX"))
	       )
	 (file-close fd1)

         (let ((data1 (sort
                       (fold 
                        (lambda (vs ax)
                          (fold
                           (lambda (v ax)
                             (let ((input-resistance (f64vector-ref v 1)))
                               (cons input-resistance ax)))
                           ax vs))
                        '() data) <))
               (dataport (open-output-file temp-path1)))
           (for-each (lambda (x) (fprintf dataport "~A~%" x)) data1)
	   (close-output-port dataport))
	 
	 (plot:init 'png (make-pathname
                          "." 
                          (sprintf "~A_results.png" 
                                   (pathname-strip-directory
                                    (pathname-strip-extension plot-label )))))
	 
	 (plot:arg "-cm" )
	 (plot:arg "-pagesize"   "12,20");;PAPER
	 (plot:arg "-textsize"   "12")
	 (plot:arg "-cpulimit"   "360")
	 (plot:arg "-maxrows"    "2001000")
	 (plot:arg "-maxfields"  "2800000")
	 (plot:arg "-maxvector"  "700000")
	 
	 (plot:proc "getdata"
		  `(
;		    ("showdata"   . "yes")
		    ("delim"      . "comma")
		    ("fieldnames" . "inputR")
		    ("pathname"   . ,temp-path1)
		    ))

	 (plot:proc "processdata"
		  `(
		    ("showdata"   . "yes")
		    ("action"      . "count")
                    ("binsize"    . 10)
		    ("fields" . "inputR")
		    ("fieldnames" . "inputRbin inputRcount")
		    ))


	 (plot:proc "areadef"
		  `(("title"     . ,(sprintf "~A" plot-label))
                    ("titledetails" . "adjust=0,0.2")
		    ("rectangle" . "2 3.5 8 10.5")
;;		    ("rectangle" . "2 5 10 14")
		    ("areacolor" . "white")

		    ("xrange"          . "0 1000")
		    ("xaxis.axisline"  . "no")
		    ("xaxis.tics"      . "no")
		    ("xaxis.stubs"     . "inc 200")
		    ("xaxis.stubrange" . "0")
;;		    ("xaxis.stubdetails" . "adjust=0,1")

		    ("yautorange"      . "datafield=inputRcount lowfix=0")
;;		    ("yaxis.label"     . "Neuron #")
		    ("yaxis.axisline"  . "no")
		    ("yaxis.tics"      . "no")
		    ("yaxis.stubs"     . "inc 10000")
;;		    ("yaxis.stubrange" . "0")
		    )
		  )
		    
       (plot:proc "bars"
		  `(("locfield"    .  "inputRbin")
		    ("lenfield"    .  "inputRcount")
		    ("thinbarline"    .  "color=gray(0.5)")
                    ))
       
       (plot:end)

       ))
)

(apply tbl-plot (command-line-arguments))
