
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
    (number-of-spikes . 16)
    (mean-firing-rate . 17)
    (mean-isi . 19)
    (isi-adaptation4 . 23)
    )
  )


(define (select-data data varname #!key (filt identity))
  (let ((var-index (alist-ref varname data-indices)))
    (match-let (((xlst xmin xmax xsum n)
                 (fold 
                  (lambda (vs ax)
                    (fold
                     (match-lambda* ((v (lst xmin xmax xsum n))
                                     (let ((x (f64vector-ref v var-index)))
                                       (if (filt x)
                                           (list (cons x lst)
                                                 (min x xmin)
                                                 (max x xmax)
                                                 (+ x xsum)
                                                 (+ n 1))
                                           (list lst
                                                 xmin
                                                 xmax
                                                 xsum
                                                 n)
                                           ))))
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


(define (output-data-range data data-path #!key (groups #f ))
  (let ((dataport (open-output-file data-path)))
    (let recur ((groups (or groups (list-tabulate (length data) (lambda (i) (+ i 1))))) 
                (data data))
      (if (null? data)
          (close-output-port dataport)
          (let ((group (car groups)))
            (for-each (lambda (x) (fprintf dataport "~A,~A~%" group x)) (car data))
            (recur (cdr groups) (cdr data)))
          ))
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
                                        
               ("titledetails" . "adjust=-0,0.2")
               ("rectangle" . ,(sprintf "~A ~A ~A ~A" x y (+ 12 x) (+ 10 y)))
               ("areacolor" . "white")

               ("xrange"          . ,(sprintf "~A ~A" xmin xmax))
               ("xaxis.axisline"  . "no")
               ("xaxis.tics"      . "yes")
               ("xaxis.stubs"     . ,(sprintf "inc ~A" xinc))
               ;("xaxis.stubrange" . "0")
  	       ("xaxis.label"     . ,xlabel)
  	       ("xaxis.labeldetails" . "adjust=-0.1,-1.3")
               ;"xaxis.stubdetails" . "adjust=0,-0.1")

               ("yautorange"      . "datafield=inputCount lowfix=0")
  	       ("yaxis.label"     . ,ylabel)
               ("yaxis.axisline"  . "no")
               ("yaxis.tics"      . "yes")
               ("yaxis.stubs"     . ,(sprintf "inc ~A" yinc))
               ("yaxis.stubrange" . "0")
               ("yaxis.stubdetails"  . "adjust=-0.1,0,size=16")
  	       ("yaxis.labeldetails" . "adjust=-1.0,0")
               )
             )
		    
  (plot:proc "bars"
             `(("locfield"    .  "inputBin")
               ("lenfield"    .  "inputCount")
               ("color"       .  "oceanblue")
               ;("thinbarline"    .  "color=oceanblue")
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
                                        
               ("titledetails" . "adjust=0,1.8")
               ("rectangle" . ,(sprintf "~A ~A ~A ~A" x y (+ 12 x) (+ 10 y)))
               ("areacolor" . "white")

               ("xautorange"      . "datafield=1 lowfix=0")
               ("xaxis.axisline"  . "no")
               ("xaxis.tics"      . "yes")
               ("xaxis.stubs"     . ,(sprintf "inc ~A" xinc))
               ("xaxis.stubrange" . "0")
  	       ("xaxis.label"     . ,xlabel)
  	       ("xaxis.labeldetails" . "adjust=-0.2,-1.3")
               ;"xaxis.stubdetails" . "adjust=0,-0.1")

               ("yautorange"      . "datafield=2")
  	       ("yaxis.label"     . ,ylabel)
               ("yaxis.axisline"  . "no")
               ("yaxis.tics"      . "yes")
               ("yaxis.stubs"     . ,(sprintf "inc ~A" yinc))
               ("yaxis.stubdetails"  . "adjust=-0.1,0,size=16")
  	       ("yaxis.labeldetails" . "adjust=-1.6,0")
               )
             )



  (plot:proc "processdata"
             `(
               ("showdata"    . "yes")
               ("fields"      . "1")
               ("action"      . "summary")
               ("valfield"    . 2)
               ))
  
  
  (plot:proc "boxplot"
             `(
               ("locfield"    . 1)
               ("basis"       . "mean")
               ("color"       . "oceanblue")
               ("printn"      , "no")
               ("mediansym"   . "shape=circle style=fill fillcolor=pink radius=0.03")
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
                     ((udata11 xmin11 xmax11 xmean11) (select-data data 'number-of-spikes))
                     ((udata12 xmin12 xmax12 xmean12) (select-data data 'mean-firing-rate))
                     ((udata15 xmin15 xmax15 xmean15) (select-data data 'isi-adaptation4
                                                                   filt: (lambda (x) (<= x 1.0))))
                     )

                    (let (
                          (data1 (sort udata1 <))
                          (data2 (sort udata2 <))
                          (data3 (sort udata3 <))
                          (data4 (sort udata4 <))
                          (data5 (sort udata5 <))
                          (data11 (sort udata11 <))
                          (data15 (sort udata15 <))
                          )

                    (plot:init 'eps (make-pathname
                                     "." 
                                     (sprintf "~A_results.eps" 
                                              (pathname-strip-directory
                                               (pathname-strip-extension plot-label )))))
                    
                    (plot:arg "-cm" )
                    (plot:arg "-textsize"   "18")
                    (plot:arg "-cpulimit"   "600")
                    (plot:arg "-maxrows"    "5001000")
                    (plot:arg "-maxfields"  "20000000")
                    (plot:arg "-maxvector"  "700000")

                    (output-data data1 temp-path1)
                    (plot-hist 3 3.5 temp-path1 1
                               "Input Resistance" "Input Resistance [MOhm]" 
                               200 0 xmax1 xmean1 "" 5000)
                    (plot:proc "page" '())

                    (output-data data2 temp-path1)
                    (plot-hist 3 3.5 temp-path1 1
                               "Time Constant" "Time Constant [ms]" 
                               5 0 xmax2 xmean2 "" 50000)
                    (plot:proc "page" '())

                    (output-data data3 temp-path1)
                    (plot-hist 3 3.5 temp-path1 1
                               "AP Threshold" "AP Threshold [mV]" 
                               4 (floor xmin3) (ceiling xmax3) xmean3 "" 50000)
                    (plot:proc "page" '())

                    (output-data data4 temp-path1)
                    (plot-hist 3 3.5 temp-path1 1
                               "AP Amplitude" "AP Amplitude [mV]" 
                               20 (floor xmin4) (ceiling xmax4) xmean4 "" 10000)
                    (plot:proc "page" '())

                    (output-data data5 temp-path1)
                    (plot-hist 3 3.5 temp-path1 1
                               "Fast AHP" "Fast AHP [mV]" 
                               10 (floor xmin5) (ceiling xmax5) xmean5 "" 50000)
                    (plot:proc "page" '())


                    (output-data-range (list udata6 udata7 udata8 udata9 udata10) 
                                       temp-path1 
                                       groups: (list 50 100 150 200 250) )
                    (plot-rangebar 3 3.5 temp-path1 2
                                   "Dendritic Attenuation" 
                                   "Distance from Soma [um]" 50
                                   "Relative AP Amplitude" 0.2)
                    (plot:proc "page" '())

                    
                    ;(plot-hist x y data-path binsize plot-label xlabel xinc xmin xmax xmean ylabel yinc)

                    (output-data data11 temp-path1)
                    (plot-hist 3 3.5 temp-path1 2
                               "AP Count" "# AP" 
                               20 (floor xmin11) xmax11 xmean11 "" 10000)
                    (plot:proc "page" '())

                    (output-data data15 temp-path1)
                    (plot-hist 2 3.5 temp-path1 0.05
                               "ISI Adaptation" "Adaptation b/w first and last AP" 
                               0.2 (floor xmin15) (min (+ xmax15 (* 0.1 xmax15)) 1.0) xmean15 "" 50000)

                    (plot:end)
         
                    ))
  ))
  )

(apply tbl-plot (command-line-arguments))
