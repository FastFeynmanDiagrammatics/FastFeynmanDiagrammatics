namespace ffd::user_space{

  bool
  operator<(QF const& q1,
	    QF const& q2){
    return ((q1.component !=
	     q2.component) ?
	    (q1.component <
	     q2.component) :
	    ((q1.statistics !=
	      q2.statistics) ?
	     (q1.statistics <
	      q2.statistics):
	     ((q1.direction !=
	       q2.direction) ?
	      (q1.direction <
	       q2.direction) :
	      (q1.is_nambu <
	       q2.is_nambu))));
  }

}//namespace
