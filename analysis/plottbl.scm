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
           (inexact->exact (truncate (round (* 100.0 (abs (- x (truncate x)))))))))

(define (sample n v)
  (let ((ub (vector-length v)))
    (let ((idxs (list-tabulate n (lambda (i) (random ub)))))
      (map (lambda (i) (vector-ref v i)) idxs)
      ))
    )


(define data-indices
  `(
    (input-resistance . 1)
    (membrane-tau . 5)
    (ap-threshold . 8)
    (ap-amplitude . 9)
    (ap-ahp . 10)
    (rel-amplitude-dend1 . 11)
    (rel-amplitude-dend2 . 12)
    (rel-amplitude-dend3 . 13)
    (rel-amplitude-dend4 . 14)
    (rel-amplitude-dend5 . 15)
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
               (list xlst xmin xmax (/ xsum n))
               ))
    )



(define (output-data data data-path)
  (let ((dataport (open-output-file data-path)))
    (for-each (lambda (x) (fprintf dataport "~A~%" x)) data)
    (close-output-port dataport)
    ))


(define (output-data-range data data-path)
  (let ((n (length (car data)))
        (dataport (open-output-file data-path)))
    (let recur ((i 0) (data data))
      (if (< i n)
          (let ((row (map car data)))
            (for-each (lambda (x) (fprintf dataport "~A" x)) (intersperse row ","))
            (fprintf dataport "~%")
            (recur (+ i 1) (map cdr data)))
          ))
    (close-output-port dataport)
    ))



(define (plot-hist x y data-path binsize plot-label xlabel xinc xmin xmax xmean ylabel yinc)
	 
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
               ("binsize"     . ,binsize)
               ("fields"      . "input")
               ("fieldnames"  . "inputBin inputCount")
               ))
  

  (plot:proc "areadef"
             `(("title"     . ,(sprintf "~A\n(Mean: ~A)" 
                                        plot-label (decimal-string xmean)))
                                        
               ("titledetails" . "adjust=0,0.2")
               ("rectangle" . ,(sprintf "~A ~A ~A ~A" x y (+ 6 x) (+ 7 y)))
               ("areacolor" . "white")

               ("xrange"          . ,(sprintf "~A ~A" xmin xmax))
               ("xaxis.axisline"  . "no")
               ("xaxis.tics"      . "yes")
               ("xaxis.stubs"     . ,(sprintf "inc ~A" xinc))
               ;("xaxis.stubrange" . "0")
  	       ("xaxis.label"     . ,xlabel)
  	       ("xaxis.labeldetails" . "adjust=0,-0.3")
               ;"xaxis.stubdetails" . "adjust=0,-0.1")

               ("yautorange"      . "datafield=inputCount lowfix=0")
  	       ("yaxis.label"     . ,ylabel)
               ("yaxis.axisline"  . "no")
               ("yaxis.tics"      . "yes")
               ("yaxis.stubs"     . ,(sprintf "inc ~A" yinc))
               ("yaxis.stubrange" . "0")
               ("yaxis.stubdetails"  . "adjust=-0.1,0")
  	       ("yaxis.labeldetails" . "adjust=-0.5")
               )
             )
		    
  (plot:proc "bars"
             `(("locfield"    .  "inputBin")
               ("lenfield"    .  "inputCount")
               ("color"       .  "oceanblue")
               ;("thinbarline"    .  "color=gray(0.5)")
               ))
  )


(define (plot-rangebar x y data-path n plot-label xlabel xinc ylabel yinc)
	 
  (plot:proc "getdata"
             `(
               ("delim"      . "comma")
               ("pathname"   . ,data-path)
               ))

  (plot:proc "areadef"
             `(("title"     . ,(sprintf "~A" plot-label))
                                        
               ("titledetails" . "adjust=0,0.2")
               ("rectangle" . ,(sprintf "~A ~A ~A ~A" x y (+ 6 x) (+ 7 y)))
               ("areacolor" . "white")

               ("xrange"          . "0 2")
               ("xaxis.axisline"  . "no")
               ("xaxis.tics"      . "yes")
               ;("xaxis.stubs"     . ,(sprintf "inc ~A" xinc))
               ("xaxis.stubrange" . "0")
  	       ("xaxis.label"     . ,xlabel)
  	       ("xaxis.labeldetails" . "adjust=0,-0.3")
               ;"xaxis.stubdetails" . "adjust=0,-0.1")

               ("yautorange"      . "datafield=1,2")
  	       ("yaxis.label"     . ,ylabel)
               ("yaxis.axisline"  . "no")
               ("yaxis.tics"      . "yes")
               ;("yaxis.stubs"     . ,(sprintf "inc ~A" yinc))
               ("yaxis.stubdetails"  . "adjust=-0.1,0")
  	       ("yaxis.labeldetails" . "adjust=-0.5")
               )
             )


  (plot:proc "categories" 
             '(("axis" . "x") 
               ("slideamount" . "-2")))


  (let recur ((i 1))
    (if (<= i n)
        (begin
	 
          (plot:proc "processdata"
                     `(
                       ("showdata"    . "yes")
                       ("action"      . "summary")
                       ("valfield"    . ,i)
                       ))
  

          (plot:proc "boxplot"
                     `(
                       ("basis"       .  "mean")
                       ("color"       .  "oceanblue")
                       ("printn"      , "no")
                       ("mediansym"   . "shape=circle style=fill fillcolor=pink radius=0.03")
                     ))

          (plot:proc "usedata" '())

          (recur (+ i 1)))
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

         (match-let (
                     ((udata1 xmin1 xmax1 xmean1) (select-data data 'input-resistance))
                     ((udata2 xmin2 xmax2 xmean2) (select-data data 'membrane-tau))
                     ((udata3 xmin3 xmax3 xmean3) (select-data data 'ap-threshold))
                     ((udata4 xmin4 xmax4 xmean4) (select-data data 'ap-amplitude))
                     ((udata5 xmin5 xmax5 xmean5) (select-data data 'ap-ahp))
                     ((udata6 xmin6 xmax6 xmean6) (select-data data 'rel-amplitude-dend1))
                     ((udata7 xmin7 xmax7 xmean7) (select-data data 'rel-amplitude-dend2))
                     ((udata8 xmin8 xmax8 xmean8) (select-data data 'rel-amplitude-dend3))
                     ((udata9 xmin9 xmax9 xmean9) (select-data data 'rel-amplitude-dend4))
                     ((udata10 xmin10 xmax10 xmean10) (select-data data 'rel-amplitude-dend5))
                     )

                    (let (
                          (data1 (sort udata1 <))
                          (data2 (sort udata2 <))
                          (data3 (sort udata3 <))
                          (data4 (sort udata4 <))
                          (data5 (sort udata5 <))
                          )

                    (plot:init 'eps (make-pathname
                                     "." 
                                     (sprintf "~A_results.eps" 
                                              (pathname-strip-directory
                                               (pathname-strip-extension plot-label )))))
                    
                    (plot:arg "-cm" )
                    (plot:arg "-pagesize"   "21.5,27.9");;PAPER
                    (plot:arg "-textsize"   "12")
                    (plot:arg "-cpulimit"   "500")
                    (plot:arg "-maxrows"    "2001000")
                    (plot:arg "-maxfields"  "3000000")
                    (plot:arg "-maxvector"  "700000")

                    (output-data data1 temp-path1)
                    (plot-hist 2 3.5 temp-path1 10
                               "Input Resistance" "Input Resistance [MOhm]" 
                               100 0 xmax1 xmean1 "Cell count" 10000)

                    (output-data data2 temp-path1)
                    (plot-hist 10 3.5 temp-path1 1 
                               "Membrane Time Constant" "Membrane Time Constant [ms]" 
                               5 0 xmax2 xmean2 "" 50000)

                    (output-data data3 temp-path1)
                    (plot-hist 2 24 temp-path1 1 
                               "AP Threshold" "AP Threshold [mV]" 
                               4 (floor xmin3) (ceiling xmax3) xmean3 "" 50000)

                    (output-data data4 temp-path1)
                    (plot-hist 10 24 temp-path1 1 
                               "AP Amplitude" "AP Amplitude [mV]" 
                               20 (floor xmin4) (ceiling xmax4) xmean4 "" 10000)

                    (output-data data5 temp-path1)
                    (plot-hist 4 14 temp-path1 1 
                               "Fast AHP" "Fast AHP [mV]" 
                               5 (floor xmin5) (ceiling xmax5) xmean5 "" 50000)

                    
                    (plot:proc "page" '())

                    (output-data-range (list udata6 udata7) temp-path1)
                    (plot-rangebar 2 3.5 temp-path1 2
                               "Dendritic Attenuation" 
                               "Dendritic AP Amplitude [mV]" 10
                               "Cell count" 10000)
         
                    (plot:end)
         
                    ))
  ))
  )

(apply tbl-plot (command-line-arguments))
