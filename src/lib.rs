//#[deny(bad_style,
//       dead-code,
//       improper-ctypes,
//       non-shorthand-field-patterns,
//       no-mangle-generic-items,
//       overflowing-literals,
//       path-statements ,
//       patterns-in-fns-without-body,
//       private-in-public,
//       unconditional-recursion,
//       unused,
//       unused-allocation,
//       unused-comparisons,
//       unused-parens,
//       while-true)]
#![allow(non_snake_case)]
#[allow(unused)]
#[allow(dead_code)]


//pub mod astronomicals
//{
	const THRESHOLD: i32 = 10;
	const MARCH: i32 = 1;
	const JUNE:i32 = 2;
	const SEPTEMBER:i32 = 3;
	const DECEMBER:i32 = 4;
	const PI:f64 = 3.1415926535897932384626433832795028841971;
	const DEGREESTORADIANS:f64 = 180.0 / PI;
	const FIRST_GREGORIAN_YEAR:i32 = 1583;
	const SECONDS_PER_DAY:f64 = 86400.0;
	const MINUTES_PER_DAY:f64 = 1440.0;
	const HOURS_PER_DAY:f64 = 24.0;



	pub struct Ephermeris
	{
		pub march_equinox: bool,// = false;
		pub june_solstice: bool,//i32 = 0;
		pub september_equinox: bool,// i32 = 0;
		pub december_solstice: bool,// i32 = 0;
		pub solstice: bool,// i32 = 0;
		pub equinox: bool,// i32 = 0;
		pub easter:bool,
	}

	fn floor(f:f64) -> f64
	{
		return (f as i64) as f64
	}

	fn dtor(degrees:f64) ->f64
	{
		return degrees / DEGREESTORADIANS
	}	


	fn cosd(d: f64) -> f64
	{
		return dtor(d).cos()
	}

	fn calcualteS(T: f64) -> f64
	{
		return 
		485.0  *cosd(324.96	+  ( 1934.136	*T))+
		203.0  *cosd(337.23	+  (32964.467	*T))+
		199.0  *cosd(342.08	+  (   20.186	*T))+
		182.0  *cosd( 27.85	+ (445267.112	*T))+
		156.0  *cosd( 73.14	+  (45036.886	*T))+
		136.0  *cosd(171.52	+  (22518.443	*T))+
		 77.0  *cosd(222.54	+  (65928.934	*T))+
		 74.0  *cosd(296.72	+  ( 3034.906	*T))+
		 70.0  *cosd(243.58	+  ( 9037.513	*T))+
		 58.0  *cosd(119.81	+  (33718.147	*T))+
		 52.0  *cosd(297.17	+  (  150.678	*T))+
		 50.0  *cosd( 21.02	+  ( 2281.226	*T))+
		 45.0  *cosd(247.54	+  (29929.562	*T))+
		 44.0  *cosd(325.15	+  (31555.956	*T))+
		 29.0  *cosd( 60.93	+  ( 4443.417	*T))+
		 18.0  *cosd(155.12	+  (67555.328	*T))+
		 17.0  *cosd(288.79	+  ( 4562.452	*T))+
		 16.0  *cosd(198.04	+  (62894.029	*T))+
		 14.0  *cosd(199.76	+  (31436.921	*T))+
		 12.0  *cosd( 95.39	+  (14577.848	*T))+
		 12.0  *cosd(287.11	+  (31931.756	*T))+
		 12.0  *cosd(320.81	+  (34777.259	*T))+
		  9.0  *cosd(227.73	+  ( 1222.114	*T))+
		  8.0  *cosd( 15.45	+  (16859.074	*T))
		
	}


	fn calculateJDE(Y_param:f64, mode: i32, debug:bool) -> f64
	{
		let Y:f64 = (Y_param-2000.0)/1000.0;
		let Y2: f64 = Y*Y;
		let Y3:f64 = Y2*Y;
		let Y4: f64 = Y3*Y;
		let mut JDEo: f64 = 0.0;
		let mut sMode = "UNKNOWN";

		match mode
		{
			MARCH =>
			{
				JDEo = 2451623.80984 + 365242.37404 * Y + 0.05169*Y2 - 0.00411*Y3 - 0.00057*Y4;
				sMode = "MARCH";
			}
			JUNE =>
			{
				JDEo = 2451716.56767 + 365241.62603 * Y + 0.00325*Y2 + 0.00888*Y3 - 0.00030*Y4;
				sMode = "JUNE";
			}
			SEPTEMBER =>
			{
				JDEo = 2451810.21715 + 365242.01767 * Y - 0.11575*Y2 + 0.00337*Y3 + 0.00078*Y4;
				sMode = "SEPTEMBER";
			}
			DECEMBER =>
			{
				JDEo = 2451900.05952 + 365242.74049 * Y - 0.06223*Y2 - 0.00823*Y3 + 0.00032*Y4;
				sMode = "DECEMBER";
			}
			_ =>
			{
				println!("Unknown mode: \"{}\"", mode);
			}
		};

		let T:f64 = (JDEo-2451545.0) / 36525.0;
		let W:f64 = (35999.373*T)-2.47;
		let DL:f64 = 1.0+(0.0334*cosd(W)) + (0.0007*cosd(2.0*W));
		let S:f64 = calcualteS(T);
		let JDE:f64 = JDEo + (0.00001*S) / DL;

		if debug
		{
			println!("\nMode: {}", sMode);
			println!("Y:      {}", Y);
			println!("JDEo:   {}", JDEo);
			println!("T:      {}", T);
			println!("W:      {}", W);
			println!("DL:     {}", DL);
			println!("S:      {}", S);
			println!("JDE:    {}", JDE);
			printDateFromJD(JDE);
		}
		return JDE
	}
	//
	//	calculateJD	calculates JD based on the given date. D, day, can include a fractional part.
	//	calculateJD assumes all dates are Gregorian. 
	//
	pub fn calculateJD(mut Y:f64, mut M:f64, D:f64) ->f64
	{
		Y = if M<3.0 {Y-1.0} else {Y};
		M = if M<3.0 {M+12.0} else{M};
		let A:f64 = floor(Y/100.0);
		let B:f64 = 2.0-A+floor(A/4.0);
		//printf("D: %lf\n", D);
		// calculate julian day for the date passed in
		let JD:f64 = floor(365.25*(Y+4716.0)) + floor(30.6001 * (M+1.0)) + D + B - 1524.5;
		//alog("\ncalculateJD\nY:%lf\n", Y);
		//alog("M:%lf\n", M);
		//alog("D:%lf\n", D);
		//alog("A:%lf\n", A);
		//alog("B:%lf\n", B);
		//alog("JD:%lf\n\n", JD);
		return JD
	}

	pub fn calculateJDLong(year:i32, month:i32, day:i32, hour:i32, minute:i32, second:i32) -> f64
	{
		let	Y:f64 =year as f64;
		let	M:f64 =month as f64;
		let	D:f64 =day as f64 + (hour as f64 / HOURS_PER_DAY) 
								+ (minute as f64 / MINUTES_PER_DAY) 
								+ (second as f64 / SECONDS_PER_DAY);
		return calculateJD(Y,M,D)
	}

	pub fn printDateFromJD(JD:f64)
	{
		let Z:f64 = ((JD +0.5) as i64) as f64;
		let mut F:f64 = (JD+0.5)-Z;
		let mut A:f64 = Z;
		if Z>=2299161.0
		{
			let alpha:f64 = floor((Z-1867216.25) / 36524.25);
			A = Z + 1.0 + (alpha - floor(alpha/4.0));
		}

		let B:f64 = A + 1524.0;
		let C:f64 = floor((B-122.1)/365.25);
		let D:f64 = floor(365.25*C);
		let E:f64 = floor((B-D) / 30.6001);

		let JD_DOM:f64 = B - D - floor(30.6001 * E) + F;
		let JD_MOY:f64 = if E< 14.0 {E-1.0}else{E - 13.0};
		let JD_YEAR:f64 = if JD_MOY>2.0 {C-4716.0}else{C-4715.0};

		let H:f64 = ((F * 24.0) as i64) as f64;
		F = F - (H / 24.0);
		let M:f64 = ((F * 24.0*60.0) as i64) as f64;
		F = F - (M / (24.0*60.0));
		let S:f64 = ((F * 24.0*60.0*60.0) as i64) as f64;


		println!("JD Date: {}-{}-{} {}:{}:{}", JD_YEAR as i32, JD_MOY as i32, JD_DOM as i32, H as i32, M as i32, S as i32);
	}


	//
	// calculateEaster calculations adapted from "Astronomical Algorithms", second edition, 2009 corrections
	// by Jean Meeus, chapter 8, pp 67-69. 
	//
	//JULIANDATE calculateEaster(int year)
	pub fn calculateEaster(year:i32) -> f64
	{
		let JD:f64;
		if year >= FIRST_GREGORIAN_YEAR
		{
			let a:i32 = year % 19;
			let b:i32 = year / 100;
			let c:i32 = year % 100;
			let d:i32 = b / 4;
			let e:i32 = b % 4;
			let f:i32 = (b+8) / 25;
			let g:i32 = (b-f+1) / 3;
			let h:i32 = ((19*a)+b-d-g+15) % 30 ;
			let i:i32 = c / 4;
			let k:i32 = c % 4;
			let l:i32 = (32+(2*e)+(2*i)-h-k) % 7;
			let m:i32 = (a+(11*h)+(22*l)) / 451;
			let np:i32 = h+l-(7*m)+114;
			let n:i32 = np / 31;
			let p:i32 = np % 31;
			JD = calculateJDLong(year, n, p+1, 0, 0, 0);
		}
		else 
		{
			let a:i32 = year % 4;
			let b:i32 = year % 7;
			let c:i32 = year % 19;
			let d:i32 = ((19*c)+15) % 30;
			let e:i32 = ((2*a)+(4*b)-d+34) % 7;
			let f:i32 = (d+e+114) / 31;
			let g:i32 = (d+e+114) % 31;
			JD = calculateJDLong(year, f, g+1,0,0,0);
		}
		return JD
	}

	pub fn calculateEphermeris(year:i32, month:i32, day:i32, _doy:i32, debug:bool) -> Ephermeris
	{
		let mut march_equinox = false;
		let mut june_solstice = false;
		let mut september_equinox = false;
		let mut december_solstice = false;
		let mut solstice = false;
		let mut equinox = false;
		let mut easter = false;

		// calcualtions are approximate, and only valid for years 1000-3000CE
		let Y:f64=year as f64;
		let M:f64 = month as f64;
		let D:f64 = day as f64;
		let jdY:f64 = if M<3.0 {(year-1) as f64} else {Y};
		let jdM:f64 = if M<3.0 {M+12.0}else{M};
		let A:f64 = floor(jdY/100.0);
		let B:f64 = 2.0-A+floor(A/4.0);

		// calcualte julian day for the date passed in
		let mut JD:f64 = floor(365.25*(jdY+4716.0)) + floor(30.6001 * (jdM+1.0)) + D + B - 1524.5;

		// calculate julian days for each for the 4 epheremis events for the year.
		let mut JDEm:f64 = calculateJDE(year as f64, MARCH, debug);
		let mut JDEj:f64 = calculateJDE(year as f64, JUNE, debug);
		let mut JDEs:f64 = calculateJDE(year as f64, SEPTEMBER, debug);
		let mut JDEd:f64 = calculateJDE(year as f64, DECEMBER, debug);
		let mut JDEe:f64 = calculateEaster(year);
		
		// now calculate whether today is easter, or not.
		JDEe -= JD;
		if (JDEe >= 0.0) && (JDEe < 1.0)
		{
			easter=true;
			//println!("JD:    {}", JD);
			//println!("Easter:{}", JDEe);
			//println!("Easter:{}", easter);
		}


		// now calcualte whether or not today (JD) is one of these ephermeris days
		// this is done by adding 0.5 to the JD and subtracting that from the JDEx.
		// if the difference is >= 0.0 and < 1.0, then the epheremis is on that day.
		JD += 0.5;
		JDEm -= JD;
		JDEj -= JD;
		JDEs -= JD;
		JDEd -= JD;
		

		if (JDEm >= 0.0) && (JDEm < 1.0)
		{
			march_equinox = true;
			equinox=true;
		}
		if (JDEj >= 0.0) && (JDEj < 1.0)
		{
			june_solstice = true;
			solstice=true;
		}
		if (JDEs >= 0.0) && (JDEs < 1.0)
		{
			september_equinox = true;
			equinox=true;
		}
		if (JDEd >= 0.0) && (JDEd < 1.0)
		{
			december_solstice = true;
			solstice=true;
		}

		let ephermeris = Ephermeris{march_equinox, june_solstice, september_equinox,
										december_solstice,solstice,equinox,easter};
		//printf("JD:     %lf\n", JD);
		//printf("MARCHEQUINOX     = %d\n", MARCHEQUINOX);
		//printf("JUNESOLSTICE     = %d\n", JUNESOLSTICE);
		//printf("SEPTEMBEREQUINOX = %d\n", SEPTEMBEREQUINOX);
		//printf("DECEMBERSOLSTICE = %d\n", DECEMBERSOLSTICE);
		return ephermeris
	}
//}