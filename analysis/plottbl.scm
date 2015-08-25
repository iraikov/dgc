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

(define (decimal-string x)
  (sprintf "~A.~A"
           (inexact->exact (truncate x))
           (inexact->exact (truncate (round (* 100.0 (- x (truncate x))))))))

(define (sample n v)
  (let ((ub (vector-length v)))
    (let ((idxs (list-tabulate n (lambda (i) (random ub)))))
      (map (lambda (i) (vector-ref v i)) idxs)
      ))
    )


(define data-indices
  `(
    (input-resistance . 1)
    (membrane-tau . 6)
    (spike-threshold . 9)
    (spike-amplitude . 10)
    (spike-ahp . 11)
    (rel-amplitude-dend1 . 12)
    (rel-amplitude-dend2 . 13)
    (rel-amplitude-dend3 . 14)
    (rel-amplitude-dend4 . 15)
    (rel-amplitude-dend5 . 16)
    (number-of-spikes . 17)
    (mean-firing-rate . 18)
    (mean-isi . 19)
    (isi-adaptation3 . 23)
    (isi-adaptation4 . 24)
    )
  )


(define (select-data data varname)
  (let ((var-index (alist-ref varname data-indices)))
    (match-let (((xlst xmin xmax xsum n)
                 (fold 
                  (lambda (vs ax)
                    (fold
                     (match-lambda* ((v (lst xmin xmax xsum n))
                                     (let ((x (f64vector-ref v var-index)))
                                       (list (cons x lst)
                                             (min x xmin)
                                             (max x xmax)
                                             (+ x xsum)
                                             (+ n 1)))))
                     ax vs))
                  '(() +inf.0 -inf.0 0.0 0) data)))
               (list (sort xlst <) xmin xmax (/ xsum n))
               ))
    )



(define (output-data data data-path)
  (let ((dataport (open-output-file data-path)))
    (for-each (lambda (x) (fprintf dataport "~A~%" x)) data)
    (close-output-port dataport)
    ))



(define (plot-data data-path plot-label xlabel xmin xmax xmean)
	 
  (plot:proc "getdata"
             `(
              ;("showdata"   . "yes")
               ("delim"      . "comma")
               ("fieldnames" . "input")
               ("pathname"   . ,data-path)
               ))

  (plot:proc "processdata"
             `(
               ;("showdata"   . "yes")
               ("action"      . "count")
               ("binsize"    . 10)
               ("fields" . "input")
               ("fieldnames" . "inputBin inputCount")
               ))
  

  (plot:proc "areadef"
             `(("title"     . ,(sprintf "~A (Mean: ~A)" 
                                        plot-label (decimal-string xmean)))
                                        
               ("titledetails" . "adjust=0,0.2")
               ("rectangle" . "2 3.5 8 10.5")
               ("areacolor" . "white")

               ("xrange"          . ,(sprintf "~A ~A" xmin xmax))
               ("xaxis.axisline"  . "no")
               ("xaxis.tics"      . "yes")
               ("xaxis.stubs"     . "inc 100")
               ("xaxis.stubrange" . "0")
  	       ("xaxis.label"     . ,xlabel)
               ;("xaxis.stubdetails" . "adjust=0,1")

               ("yautorange"      . "datafield=inputCount lowfix=0")
  	       ("yaxis.label"     . "Cell count")
               ("yaxis.axisline"  . "no")
               ("yaxis.tics"      . "yes")
               ("yaxis.stubs"     . "inc 10000")
               ("yaxis.stubrange" . "0")
               ("yaxis.stubdetails"  . "adjust=-0.1,0")
  	       ("yaxis.labeldetails" . "adjust=-0.3,0")
               )
             )
		    
  (plot:proc "bars"
             `(("locfield"    .  "inputBin")
               ("lenfield"    .  "inputCount")
               ("thinbarline"    .  "color=gray(0.5)")
               ))
  )


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

         (match-let (((data1 xmin1 xmax1 xmean1) (select-data data 'input-resistance)))

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
                    (plot:arg "-maxfields"  "3000000")
                    (plot:arg "-maxvector"  "700000")

                    (output-data data1 temp-path1)
                    (plot-data temp-path1 "Input Resistance" "MOhm" 0 xmax1 xmean1)
         
                    (plot:end)
         
                    ))
  ))

(apply tbl-plot (command-line-arguments))
