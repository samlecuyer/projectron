// Copyright 2015 Sam L'ecuyer. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package projectron

import "math"
import "errors"

const (
	s_pi    float64 = 3.14159265359
	two_pi  float64 = math.Pi * 2
	half_pi float64 = math.Pi / 2
	fort_pi float64 = math.Pi / 4
	d2r     float64 = math.Pi / 180
)

//pj_msfn(double sinphi, double cosphi, double es) {
//	return (cosphi / sqrt (1. - es * sinphi * sinphi));
//}

func msfn(sinphi, cosphi, es float64) float64 {
	return (cosphi / math.Sqrt(1-es*sinphi*sinphi))
}

// pj_tsfn(double phi, double sinphi, double e) {
// 	sinphi *= e;
// 	return (tan (.5 * (HALFPI - phi)) /
// 	   pow((1. - sinphi) / (1. + sinphi), .5 * e));
// }

func tsfn(phi, sinphi, e float64) float64 {
	return math.Tan(.5*(half_pi-phi)) / math.Pow((1-sinphi)/(1+sinphi), .5*e)
}

//module.exports = function(eccent, ts) {
//   var eccnth = 0.5 * eccent;
//   var con, dphi;
//   var phi = HALF_PI - 2 * Math.atan(ts);
//   for (var i = 0; i <= 15; i++) {
//     con = eccent * Math.sin(phi);
//     dphi = HALF_PI - 2 * Math.atan(ts * (Math.pow(((1 - con) / (1 + con)), eccnth))) - phi;
//     phi += dphi;
//     if (Math.abs(dphi) <= 0.0000000001) {
//       return phi;
//     }
//   }
//   //console.log("phi2z has NoConvergence");
//   return -9999;
// };

func phi2(e, ts float64) (float64, error) {
	eth := e * 0.5
	phi := half_pi - 2*math.Atan(ts)
	var con, dphi float64
	for i := 0; i <= 15; i++ {
		con = e * math.Sin(phi)
		dphi = half_pi - 2*math.Atan(ts*(math.Pow(((1-con)/(1+con)), eth))) - phi
		phi += dphi
		if math.Abs(dphi) < epsln {
			return phi, nil
		}
	}
	return math.Inf(-1), errors.New("phi2 has no convergence")
}

// #define SPI     3.14159265359
// #define TWOPI   6.2831853071795864769
// #define ONEPI   3.14159265358979323846

// double adjlon (double lon) {
//     if (fabs(lon) <= SPI) return( lon );
//     lon += ONEPI;  /* adjust to 0..2pi rad */
//     lon -= TWOPI * floor(lon / TWOPI); /* remove integral # of 'revolutions'*/
//     lon -= ONEPI;  /* adjust back to -pi..pi rad */
//     return( lon );
// }

func adjLng(x float64) float64 {
	if math.Abs(x) <= s_pi {
		return x
	} else {
		x += math.Pi
		x -= two_pi * math.Floor(x/two_pi)
		x -= math.Pi
		return x
	}
}

func sign(x float64) float64 {
	if math.Signbit(x) {
		return -1
	} else {
		return 1
	}
}
