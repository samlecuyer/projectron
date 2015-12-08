// Copyright 2015 Sam L'ecuyer. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package projectron

import (
	"math"
	"testing"
	// "fmt"
)

func TestDegreeString(t *testing.T) {
	for _, pm := range pm_list {
		parseDegreeString(pm.defn)
	}
}

func TestProjString(t *testing.T) {
	pj, err := NewProjection("+title=WGS 84 (long/lat) +proj=longlat +ellps=WGS84 +datum=WGS84 +units=degrees")
	if err != nil {
		t.Error(err)
	}
	lng0, lat0 := 18.5*d2r, 54.2*d2r
	x, y, err := pj.Forward(lng0, lat0)
	if err != nil {
		t.Error(err)
	}

	// should translate forward
	if !close(lng0, x) || !close(lat0, y) {
		t.Errorf("fwd translation off: (%f, %f) - (%f, %f)", lng0, lat0, x, y)
	}

	lng1, lat1, err := pj.Inverse(x, y)
	if err != nil {
		t.Error(err)
	}
	// should translate back
	if !close(lng0, lng1) || !close(lat0, lat1) {
		t.Errorf("inv translation off: (%f, %f) - (%f, %f)", lng0, lat0, lng1, lat1)
	}
}

func TestMercator(t *testing.T) {
	pj, err := NewProjection("+title=WGS 84 / Pseudo-Mercator +proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs")
	if err != nil {
		t.Error(err)
	}
	lng0, lat0 := 18.5*d2r, 54.2*d2r
	expx, expy := 2059410.57968, 7208125.2609
	x, y, err := pj.Forward(lng0, lat0)
	if err != nil {
		t.Error(err)
	}

	// should translate forward
	if !close(expx, x) || !close(expy, y) {
		t.Errorf("fwd translation off: (%f, %f) - (%f, %f)", expx, expy, x, y)
	}

	lng1, lat1, err := pj.Inverse(x, y)
	if err != nil {
		t.Error(err)
	}
	// should translate back
	if !close(lng0, lng1) || !close(lat0, lat1) {
		t.Errorf("inv translation off: (%f, %f) - (%f, %f)", lng0, lat0, lng1, lat1)
	}
}

func TestLCC(t *testing.T) {
	t.Skip("not implemented")
	pj, err := NewProjection("+proj=lcc +lat_0=18 +lat_1=18 +lon_0=-77 +k_0=1 +k_0=1.0 +R=6378137")
	if err != nil {
		t.Error(err)
	}

	// c, n, rho0 float64
	// phi2, phi1 float64
	// ellips bool
	// lcc, _ := pj.(*LCC)

	lng0, lat0 := -.1396263, .4712389
	// println(lng0, lat0)
	expx, expy := 8701763.068335464, -139008.773062367
	x, y, err := pj.Forward(lng0, lat0)
	if err != nil {
		t.Error(err)
	}

	// should translate forward
	if !close(expx, x) || !close(expy, y) {
		t.Errorf("fwd translation off: (%f, %f) - (%f, %f)", expx, expy, x, y)
	}

	// lng1, lat1, err := pj.Inverse(x, y)
	// if err != nil {
	// 	t.Error(err)
	// }
	// // should translate back
	// if !close(lng0, lng1) || !close(lat0, lat1) {
	// 	t.Errorf("inv translation off: (%f, %f) - (%f, %f)", lng0, lat0, lng1, lat1)
	// }
}


func close(a, b float64) bool {
	return math.Abs(a-b) < 1.0e-5
}
