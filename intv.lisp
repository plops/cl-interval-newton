;; play with interval arithmetic

(defclass interval ()
  ((a :type double-float :reader inf :initarg :a)
   (b :type double-float :reader sup :initarg :b)))

(defun make-interval (a &optional (b a))
  (make-instance
   'interval
   :a a :b b))

(defun make-interval* (a b)
  (if (< a b)
      (make-interval a b)
      (make-interval b a)))

(defmethod mid ((v interval))
  (with-slots (a b) v
    (/ (+ b a) 2)))

(defmethod .sqrt ((v interval))
  (with-slots (a b) v
    (assert (<= a b))
      (make-instance 'interval
		     :a (sqrt a)
		     :b (sqrt b))))
(defmethod .cos ((v interval))
  (with-slots (a b) v
    (assert (<= a b))
    (make-interval* (cos a) (cos b))))
(defmethod .sin ((v interval))
  (with-slots (a b) v
    (assert (<= a b))
    (make-interval* (sin a) (sin b))))
(defmethod .1+ ((v interval))
  (with-slots (a b) v
      (assert (<= a b))
      
      (make-instance 'interval
		     :a (1+ a)
		     :b (1+ b))))
(defmethod .+ ((v interval) (w interval))
  (with-slots (a b) v
    (with-slots ((x a) (y b)) w
     (assert (<= a b))
      (assert (<= x y))
      (make-instance 'interval
		    :a (+ a x)
		    :b (+ b y)))))

(defmethod .* ((v interval) (w interval))
  (with-slots (a b) v
    (with-slots ((x a) (y b)) w
      (assert (<= a b))
      (assert (<= x y))
      (make-interval* 
		    (* a x)
		    (* b y)))))

(defmethod ./ ((v interval) (w interval))
  (with-slots (a b) v
    (with-slots ((x a) (y b)) w
      (assert (<= a b))
      (assert (<= x y))
      (make-instance 'interval
		    :a (/ a y)
		    :b (/ b x)))))

#+nil
(defclass gen-interval-fun (name ))

(defmethod .in ((s double-float) (v interval))
  (with-slots (a b) v
    (assert (<= a b))
      (<= a s b)))

(defmethod .in ((v interval) (w interval))
  (with-slots (a b) v
      (with-slots ((x a) (y b)) w
	(and (<= x a y)
	     (<= x b y)))))

(defmethod .in ((v interval) (w interval))
  (with-slots (a b) v
    (with-slots ((x a) (y b)) w
      (make-interval (max a x) (min b y)))))


(defmethod .f ((v interval))
  (.+ (.sqrt v)  (.* (.1+ v)
		     (.cos v))
      ))

(defmethod .deriv ((v interval))
  (.+ (./ (make-interval 1.0)
	  (.* (make-interval 2.0) (.sqrt v)))
      (.cos v)
      (.* (.* (make-interval -1.0) (.1+ v)
	      )
	  (.sin v))))

(defun f (x)
  (.f
   (make-interval x)))


(defmethod criter-p ((v interval))
  (and (< (sup (.* (f (inf v))
		   (f (sup v))))
	  .0)
       (not (.in 0.0 (.deriv v)))))


(defun newton (y)
  (let ((yo y))
    (when (criter-p y)
      (loop do
	   (setf yo y
		 y (.and (.- (make-interval (.mid y))
			     (./ (f (.mid y))
				 (.deriv y)))
			 y))
	 while
	   (.in y yo)))))

(newton (make-interval 2.0 3.0))
