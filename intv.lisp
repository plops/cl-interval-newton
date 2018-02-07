;; play with interval arithmetic
;; reference: uni wuppertal xsc
(declaim (optimize (safety 3) (speed 0) (debug 3)))

(defclass interval ()
  ((a :type double-float :reader inf :initarg :a)
   (b :type double-float :reader sup :initarg :b)))

(defmethod print-object ((an-interval interval) a-stream)
  (with-slots (a b) an-interval
   (format a-stream "<interval ~f ~f>" a b)))

(defun make-interval (a &optional (b a))
  (declare (type double-float a b)
	   (values interval &optional))
  (assert (<= a b))
  (make-instance
   'interval
   :a a :b b))

(defun make-interval* (a b)
  (declare (type double-float a b)
	   (values interval &optional))
  (if (< a b)
      (make-interval a b)
      (make-interval b a)))

(defmethod mid ((v interval))
  (declare 
	   (values double-float &optional))
  (with-slots (a b) v
    (/ (+ b a) 2)))

(defmethod .mid ((v interval))
  (declare (values interval &optional))(declare (values interval &optional))
  (with-slots (a b) v
    (make-interval (/ (+ b a) 2))))

(defmethod .sqrt ((v interval))
  (declare (values interval &optional))
  (with-slots (a b) v
    (assert (<= a b))
    (assert (<= 0.0 a))
    (assert (<= 0.0 b))
    (make-instance 'interval
			  :a (sqrt a)
			  :b (sqrt b))))
(defmethod .cos ((v interval))
  (declare (values interval &optional))
  (with-slots (a b) v
    (assert (<= a b))
    (make-interval* (cos a) (cos b))))
(defmethod .sin ((v interval))
  (declare (values interval &optional))
  (with-slots (a b) v
    (assert (<= a b))
    (make-interval* (sin a) (sin b))))
(defmethod .1+ ((v interval))
  (declare (values interval &optional))
  (with-slots (a b) v
      (assert (<= a b))
      
      (make-instance 'interval
		     :a (1+ a)
		     :b (1+ b))))
(defmethod .- ((v interval) (w interval))
  (declare (values interval &optional))
  (with-slots (a b) v
    (with-slots ((x a) (y b)) w
     (assert (<= a b))
      (assert (<= x y))
      (make-interval* (- a y)
		      (- b x)))))

(defmethod .+ ((v interval) (w interval))
  (declare (values interval &optional))
  (with-slots (a b) v
    (with-slots ((x a) (y b)) w
     (assert (<= a b))
      (assert (<= x y))
      (make-instance 'interval
		    :a (+ a x)
		    :b (+ b y)))))

(defmethod .* ((v interval) (w interval))
  (declare (values interval &optional))
  (with-slots (a b) v
    (with-slots ((x a) (y b)) w
      (assert (<= a b))
      (assert (<= x y))
      (make-interval* 
		    (* a x)
		    (* b y)))))
(defmethod .* ((s double-float) (v interval))
  (declare (values interval &optional))
  (with-slots (a b) v
    (assert (<= a b))
	
	      (make-interval* 
	       (* a s)
	       (* b s))))

(defmethod ./ ((v interval) (w interval))
  (declare (values interval &optional))
  (with-slots (a b) v
    (with-slots ((x a) (y b)) w
      (assert (<= a b))
      (assert (<= x y))
      (make-interval*
		     (/ a y)
		     (/ b x)))))

#+nil
(defclass gen-interval-fun (name ))

(defmethod .in-p ((s number) (v interval))
  
  (with-slots (a b) v
    (assert (<= a b))
      (<= a s b)))

(defmethod .in-entirely-p ((v interval) (w interval))
  
  (with-slots (a b) v
      (with-slots ((x a) (y b)) w
	(assert (<= a b))
	(assert (<= x y))
      
	(and (< x a y)
	     (< x b y)))))

(defmethod ./= ((v interval) (w interval))
  
  (with-slots (a b) v
      (with-slots ((x a) (y b)) w
	(assert (<= a b))
	(assert (<= x y))
      
	(or
	 (/= x a)
	 (/= b y)))))

(defmethod .and ((v interval) (w interval))
  (declare (values interval &optional))
  (with-slots (a b) v
    (with-slots ((x a) (y b)) w
      (assert (<= a b))
      (assert (<= x y))
      
      (make-interval (max a x) (min b y)))))



;; newton iteration example


(defmethod .f ((v interval))
  (declare (values interval &optional))
  (.+ (.sqrt v)  (.* (.1+ v)
		     (.cos v))
      ))
#+nil
(.f (make-interval 2d5 2d5))

(defmethod .deriv ((v interval))
  (declare (values interval &optional))
  (.+ (.+ (./ (make-interval 1d0)
	   (.* (make-interval 2d0)
	       (.sqrt v)))
	  (.cos v))
      (.* (.* (make-interval -1d0)
	      (.1+ v))
	  (.sin v))))

(defun f (x)
  (.f
   (make-interval x)))


(defmethod criter-p ((v interval))
  (and (< (sup (.* (f (inf v))
		   (f (sup v))))
	  .0)
       (not (.in-p 0.0 (.deriv v)))))


(defun newton (y)
  (let ((yo y))
    (when (criter-p y)
      (loop do
	   (let ((mid (mid y)))
	    (setf yo y
		  y (.and (.- (make-interval mid)
			      (./ (f mid)
				  (.deriv y)))
			  y)))
	   (format t "~a~%" y)
	   
	 while
	   (./= y yo))
      y)))

#+nil
(.in-entirely-p 
 (make-interval 2d0 3d0)
 (make-interval 1d0 10d0))
#+nil
(newton (make-interval 2d0 3d0))
#+nil
(criter-p (make-interval 2d0 3d0))
#+nil
(./ (.f (.mid (make-interval 2d0 3d0)))
    (.deriv (make-interval 2d0 3d0)))


;; runge kutta example
;; Y' = F(x,Y),    Y(x0)=Y0

(deftype interval-array ()
  `(simple-array interval (*)))

(defun sa* (s a)
  (declare (type double-float s)
	   (type interval-array a)
	   (values interval-array &optional))
  (let* ((n (length a))
	 (res (make-array n :element-type interval)))
    (dotimes (i n)
      (setf (aref n i) (.* s  (aref a i))))
    res))

(defun a* (a b)
  (declare
   
	   (type interval-array a b)
	   (values interval-array &optional))
  (assert (= (length a)
	     (length b)))
  (let* ((n (length a))
	 (res (make-array n :element-type interval)))
    (dotimes (i n)
      (setf (aref n i) (.* (aref a i) (aref b i))))
    res))

(defun rk-f (x y)
  (declare (type interval-array y))
  (let ((z (make-array 3 :element-type interval)))
    (setf (aref z 0) (.* (aref y 2) (aref y 3))
	  (aref z 1) (.* -1d0 (.* (aref y 1) (aref y 3)))
	  (aref z 2) (.* -.522d0 (.* (aref y 1) (aref y 2))))
    z))

(defun rk-init ()
  (let ((y (make-array 3 :element-type interval)))
    (loop for e in '(0 1 1) and i from 0 do
	 (setf (aref y i) (make-interval (* 1d0 e))))
    ;; x h y
    (values 0d0 .1d0 y)))


(defun rk-run ()
  (destructuring-bind (x h y) (rk-init)
    (let ((k (make-array (list 4 3) :element-type interval)))
      (setf (aref k 0) (* h (rk-f x y))
	    (aref k 1) (* h (rk-f (+ x (* .5 h)) (a+ y (.* .5d0 (aref k 0)))))
	    (aref k 2) (* h (rk-f (+ x (* .5 h)) (a+ y (.* .5d0 (aref k 1)))))
	    (aref k 1) (* h (rk-f (+ x h) (a+ (aref k 2) y)))
	    y (.+ y (.* (/ 6d0) (.+ (.+ (aref k 0) (aref k 3))
				    (.* 2d0 (.+ (aref k 1) (aref k 2))))))
	    x (+ x h))
      (format t
	      "~a ~a~%" x y))))
