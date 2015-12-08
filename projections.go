// Copyright 2015 Sam L'ecuyer. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package projectron

import "math"
import "errors"

type impl interface {
	Projection
	init(paramset) error
}

func lookupImpl(pin *pj) impl {
	switch pin.proj {
	case "latlong", "longlat", "latlon", "lonlat":
		return &LngLat{pin}
	case "merc":
		return &Mercator{pin}
	case "lcc":
		return &LCC{pj: pin}
	case "eqc":
		return &Equirectangular{pj: pin}
	}
	return nil
}

type LngLat struct {
	*pj
}

func (ll *LngLat) init(params paramset) error {
	ll.x0 = 0
	ll.y0 = 0
	return nil
}

func (ll *LngLat) IsLngLat() bool {
	return true
}

func (ll *LngLat) Forward(lng, lat float64) (x, y float64, err error) {
	return ll.commonFwd(lng, lat, ll.fwd)
}

func (ll *LngLat) Inverse(x, y float64) (lng, lat float64, err error) {
	return ll.commonInv(x, y, ll.inv)
}

func (ll *LngLat) fwd(lam, phi float64) (float64, float64, error) {
	x := lam / ll.a
	y := phi / ll.a
	return x, y, nil
}

func (ll *LngLat) inv(x, y float64) (lng, lat float64, err error) {
	lat = y * ll.a
	lng = x * ll.a
	return lng, lat, nil
}

type Mercator struct {
	*pj
}

func (m *Mercator) IsLngLat() bool {
	return true
}
func (m *Mercator) ToMeter() float64 {
	return m.to_meter
}
func (m *Mercator) FromGreenwich() float64 {
	return m.from_greenwich
}

func (m *Mercator) init(params paramset) error {
	var phits float64
	var isPhits bool
	if phits, isPhits = params.degree("lat_ts"); isPhits {
		phits = math.Abs(phits)
	}
	if m.e != 0 && isPhits {
		m.k0 = msfn(math.Sin(phits), math.Cos(phits), m.es)
	} else if isPhits {
		m.k0 = math.Cos(phits)
	}
	return nil
}

func (m *Mercator) Forward(lng, lat float64) (x, y float64, err error) {
	return m.commonFwd(lng, lat, m.fwd)
}

func (m *Mercator) Inverse(x, y float64) (lng, lat float64, err error) {
	return m.commonInv(x, y, m.inv)
}

func (m *Mercator) fwd(lam, phi float64) (x float64, y float64, err error) {
	if m.es != 0 {
		x = m.k0 * lam
		y = -m.k0 * math.Log(tsfn(phi, math.Sin(phi), m.e))
	} else {
		x = m.k0 * lam
		y = m.k0 * math.Log(math.Tan(fort_pi+0.5*phi))
	}
	return x, y, nil
}

func (m *Mercator) inv(x, y float64) (lng, lat float64, err error) {
	if m.es != 0 {
		lat, err = phi2(math.Exp(-y/m.k0), m.e)
		lng = x * m.k0
	} else {
		lng = x / m.k0
	lat = half_pi - 2*math.Atan(math.Exp(-y/m.k0))
	}
	
	return lng, lat, err
}


type LCC struct {
	*pj
	c, n, rho0 float64
	phi2, phi1 float64
	ellips bool
}

func (lcc *LCC) IsLngLat() bool {
	return false
}

func (ll *LCC) init(params paramset) error {
	ll.phi1, _ = params.degree("lat_1")
	if phi2, ok := params.degree("lat_2"); ok {
		ll.phi2 = phi2
	} else {
		ll.phi2 = ll.phi1
		if _, ok := params.string("lat_0"); !ok {
			ll.phi0 = ll.phi1
		}
	}
	if math.Abs(ll.phi1 + ll.phi2) <= epsln {
		return errors.New("these can't be the same")
	}
	sinphi := math.Sin(ll.phi1)
	ll.n = sinphi
	cosphi := math.Cos(ll.phi1)
	secant := math.Abs(ll.phi1 - ll.phi2) >= epsln
	ll.ellips = ll.es != 0
	if ll.ellips {
		ll.e = math.Sqrt(ll.es)
		m1 := msfn(sinphi, cosphi, ll.es)
		ml1 := tsfn(ll.phi1, sinphi, ll.e)
		if secant {
			sinphi = math.Sin(ll.phi2)
			ll.n = math.Log(m1/msfn(sinphi, math.Cos(ll.phi2), ll.es))
			ll.n /= math.Log(ml1 / tsfn(ll.phi2, sinphi, ll.e))
		}
		ll.c = m1 * math.Pow(ml1, -ll.n) / ll.n
		if math.Abs(math.Abs(ll.phi0) - half_pi) < epsln {
			ll.rho0 = 0
		} else {
			ll.rho0 = ll.c * math.Pow(tsfn(ll.phi0, math.Sin(ll.phi0), ll.e), ll.n)
		}
	} else {
		if secant {
			ll.n = math.Log(cosphi / math.Cos(ll.phi2)) /
				  	math.Log(math.Tan(fort_pi + .5 * ll.phi2) /
				  		math.Tan(fort_pi + .5 * ll.phi1))
		}
		ll.c = cosphi * math.Pow(math.Tan(fort_pi + .5 * ll.phi1), ll.n) / ll.n
		if math.Abs(math.Abs(ll.phi0) - half_pi) < epsln {
			ll.rho0 = 0
		} else {
			ll.rho0 = ll.c * math.Pow(math.Tan(fort_pi + 0.5 * ll.phi0), -ll.n)
		}
	}
	// println(ll.phi1, ll.phi2, ll.n, ll.rho0, ll.c,)
	return nil
}

func (ll *LCC) Forward(lng, lat float64) (x, y float64, err error) {
	return ll.commonFwd(lng, lat, ll.fwd)
}

func (ll *LCC) Inverse(x, y float64) (lng, lat float64, err error) {
	return ll.commonInv(x, y, ll.inv)
}

func (ll *LCC) fwd(lam, phi float64) (x float64, y float64, err error) {
	var rho float64
	if math.Abs(math.Abs(phi)-half_pi) < epsln {
		if phi*ll.n <= 0 {
			return hugeVal, hugeVal, errors.New("something")
		}
	} else {
		if ll.ellips {
			rho = ll.c * math.Pow(tsfn(phi, math.Sin(phi), ll.e), ll.n)
		} else {
			rho = ll.c * math.Pow(math.Tan(fort_pi + 0.5 * phi), -ll.n)
		}
	}
	lam *= ll.n
	x = ll.k0 * (rho * math.Sin(lam))
	y = ll.k0 * (ll.rho0 - rho*math.Cos(lam))
	return x, y, nil
}

func (ll *LCC) inv(x, y float64) (lng, lat float64, err error) {
	panic("don't call this")
}

type Equirectangular struct {
	*pj
	phi1 float64
}

func (eqc *Equirectangular) init(params paramset) error {
	eqc.phi1, _ = params.degree("lat_1")
	return nil
}

func (eqc *Equirectangular) IsLngLat() bool {
	return false
}
func (eqc *Equirectangular) Forward(lng, lat float64) (x, y float64, err error) {
	return eqc.commonFwd(lng, lat, eqc.fwd)
}

func (eqc *Equirectangular) Inverse(x, y float64) (lng, lat float64, err error) {
	return eqc.commonInv(x, y, eqc.inv)
}

func (eqc *Equirectangular) fwd(lam, phi float64) (float64, float64, error) {
	x := lam / eqc.a
	y := phi / eqc.a
	return x, y, nil
}

func (eqc *Equirectangular) inv(x, y float64) (lng, lat float64, err error) {
	lat = y * eqc.a
	lng = x * math.Cos(eqc.phi1) * eqc.a
	return lng, lat, nil
}

