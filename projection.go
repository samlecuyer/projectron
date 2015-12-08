// Copyright 2015 Sam L'ecuyer. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package projectron

import (
	"errors"
	"math"
	"strconv"
	"strings"
)

type Projection interface {
	// Forward projects lng/lat into this.  l/l are in radians
	Forward(lng, lat float64) (x, y float64, err error)
	// Inverse projects this back to lng/lat
	Inverse(x, y float64) (lng, lat float64, err error)
	IsLngLat() bool
	ToMeter() float64
	FromGreenwich() float64
	Radius() float64
}

func NewProjection(str string) (Projection, error) {
	parms := make(paramset)
	var ok bool
	for _, part := range strings.Split(str, "+") {
		param := strings.TrimSpace(part)
		if param == "" {
			continue
		}
		key, val := keyVal(param)
		parms[key] = val
	}
	pin := &pj{axis: "enu"}
	if pin.proj, ok = parms.string("proj"); !ok {
		return nil, ErrUnsupportedProj
	}
	pin.setDatum(parms)
	pin.setEllipse(parms)

	pin.aOrig = pin.a
	pin.esOrig = pin.es

	pin.e = math.Sqrt(pin.es)
	pin.ra = 1 / pin.a
	pin.oneEs = 1 - pin.es
	pin.rOneEs = 1 / pin.oneEs

	if pin.datumType == PJD_3PARAM && (pin.datumParams == nil || (pin.datumParams[0] == 0 &&
		pin.datumParams[1] == 0 && pin.datumParams[2] == 0)) &&
		pin.a == 6378137.0 && math.Abs(pin.es-0.006694379990) < 0.000000000050 { /*WGS84/GRS80*/
		pin.datumType = PJD_WGS84
	}

	pin.geoc, _ = parms.bool("geoc")
	pin.over, _ = parms.bool("over")

	if lwc, ok := parms.degree("lon_wrap"); ok {
		pin.long_wrap_set = ok
		pin.long_wrap_center = lwc
	}

	if axis, ok := parms.string("axis"); ok {
		if len(axis) != 3 {
			return nil, ErrInvalidParam
		}
		// TODO: validate (I'm in a hurry)
		pin.axis = axis
	}

	// central meridian
	pin.lam0, _ = parms.degree("lon_0")
	// central latitude
	pin.phi0, _ = parms.degree("lat_0")

	// false easting/northing
	pin.x0, _ = parms.float("x_0")
	pin.y0, _ = parms.float("y_0")

	// general scaling
	if k0, ok := parms.float("k_0"); ok {
		pin.k0 = k0
	} else if k0, ok := parms.float("k"); ok {
		pin.k0 = k0
	} else {
		pin.k0 = 1.
	}
	if pin.k0 <= 0 {
		return nil, ErrInvalidParam
	}

	// units
	var s string
	if name, ok := parms.string("units"); ok {
		if unit, ok := units_list[name]; ok {
			pin.to_meter = unit.to_meter
			pin.fr_meter = 1 / unit.to_meter
		}
	} else {
		s, _ = parms.string("to_meter")
	}
	if s != "" {
		// do something about the to_meter thing
	} else {
		pin.to_meter = 1
		pin.fr_meter = 1
	}
	// vertical units
	s = ""
	if name, ok := parms.string("vunits"); ok {
		if unit, ok := units_list[name]; ok {
			pin.to_meter = unit.to_meter
			pin.fr_meter = unit.to_meter
		}
	} else {
		s, _ = parms.string("vto_meter")
	}
	if s != "" {
		// do something about the to_meter thing
	} else {
		pin.vto_meter = pin.to_meter
		pin.vfr_meter = pin.fr_meter
	}

	// prime meridian
	if name, ok := parms.string("pm"); ok {
		if pm, ok := pm_list[name]; ok {
			pin.from_greenwich = parseDegreeString(pm.defn)
		}
	}

	imp := lookupImpl(pin)
	if imp != nil {
		imp.init(parms)
		return imp, nil
	}
	return nil, ErrUnsupportedProj
}

type pj struct {
	proj                 string
	axis                 string
	datumType            datumType
	datumParams          []float64
	catalogName          string
	a, es, e, ra         float64
	lam0, phi0, k0       float64
	x0, y0               float64
	oneEs, rOneEs        float64
	aOrig, esOrig        float64
	geoc, over           bool
	long_wrap_set        bool
	long_wrap_center     float64
	to_meter, fr_meter   float64
	vto_meter, vfr_meter float64
	from_greenwich       float64

}

func (p *pj) setDatum(params paramset) error {
	if name, ok := params.string("datum"); ok {
		if datum, ok := datums_list[name]; ok {
			params["ellps"] = datum.ellipse
			key, val := keyVal(datum.definition)
			params[key] = val
		}
	}

	if _, ok := params.string("nadgrids"); ok {
		// BUG(slecuyer): no support for nadgrids
		p.datumType = PJD_GRIDSHIFT
	} else if catalog, ok := params.string("catalog"); ok {
		p.datumType = PJD_GRIDSHIFT
		p.catalogName = catalog
		// BUG(slecuyer): we don't do anything with the catalog date
	} else if towgs84, ok := params.string("towgs84"); ok {
		parts := strings.Split(towgs84, ",")
		p.datumParams = make([]float64, 7)
		for i, f := range parts {
			p.datumParams[i], _ = strconv.ParseFloat(f, 64)
		}
		if len(parts) == 7 {
			p.datumType = PJD_7PARAM
			p.datumParams[3] *= sec2rad
			p.datumParams[4] *= sec2rad
			p.datumParams[5] *= sec2rad
			p.datumParams[6] = (p.datumParams[6] / 1000000.0) + 1
		} else {
			p.datumType = PJD_3PARAM
		}
	}
	return nil
}

func (p *pj) setEllipse(params paramset) error {
	var b float64
	if r, ok := params.float("R"); ok {
		p.a = r
	} else {
		if name, ok := params.string("ellps"); ok {
			if ellps, ok := ellipse_list[name]; ok {
				key, val := keyVal(ellps.major)
				if _, ok = params[key]; !ok {
					params[key] = val
				}
				key, val = keyVal(ellps.ell)
				if _, ok = params[key]; !ok {
					params[key] = val
				}
			}
		}
		p.a, _ = params.float("a")
		if es, ok := params.float("es"); ok {
			p.es = es
		} else if e, ok := params.float("e"); ok {
			p.es = e * e
		} else if rf, ok := params.float("rf"); ok {
			p.es = 1 / rf
			p.es = p.es * (2 - p.es)
		} else if f, ok := params.float("f"); ok {
			p.es = f * (2. - f)
		} else if b, ok = params.float("b"); ok {
			p.es = 1. - (b*b)/(p.a*p.a)
		}
		if b == 0 {
			b = p.a * math.Sqrt(1.-p.es)
		}
		if ra, ok := params.bool("R_A"); ra && ok {
			p.a *= 1. - p.es*(sixth+p.es*(ra4+p.es*ra6))
			p.es = 0.
		} else if rv, ok := params.bool("R_V"); rv && ok {
			p.a *= 1. - p.es*(sixth+p.es*(rv4+p.es*rv6))
			p.es = 0.
		} else if rg, ok := params.bool("R_g"); rg && ok {
			p.a = math.Sqrt(p.a * b)
			p.es = 0.
		} else if rh, ok := params.bool("R_h"); rh && ok {
			p.a = 2 * p.a * b / (p.a + b)
			p.es = 0.
			// BUG(slecuyer): no support for R_lat_a or R_lat_g
		}
	}
	return nil
}

func (p *pj) commonFwd(lam, phi float64, tr translator) (x, y float64, err error) {
	// println(p.to_meter, p.fr_meter, p.a)
	t := math.Abs(phi) - half_pi
	if t > epsln || math.Abs(lam) > 10 {
		return hugeVal, hugeVal, errors.New("this is way out of bounds")
	}
	if math.Abs(t) <= epsln {
		phi = math.Copysign(half_pi, phi)
	} else if p.geoc {
		phi = math.Atan(p.rOneEs * math.Tan(phi))
	}
	lam -= p.lam0
	if !p.over {
		lam = adjLng(lam)
	}
	x, y, err = tr(lam, phi)
	if err != nil {
		return hugeVal, hugeVal, err
	}
	x = p.fr_meter * (p.a*x + p.x0)
	y = p.fr_meter * (p.a*y + p.y0)
	return
}

func (p *pj) commonInv(x, y float64, tr translator) (lam, phi float64, err error) {
	if x == hugeVal || y == hugeVal {
		return hugeVal, hugeVal, errors.New("this is way out of bounds")
	}
	x = (x*p.to_meter - p.x0) * p.ra
	y = (y*p.to_meter - p.y0) * p.ra
	lam, phi, err = tr(x, y)
	if err != nil {
		return hugeVal, hugeVal, err
	}
	y += p.lam0
	if !p.over {
		x = adjLng(x)
	}
	if p.geoc && math.Abs(math.Abs(phi)-half_pi) > epsln {
		phi = math.Atan(p.oneEs * math.Tan(phi))
	}
	return
}

func (p *pj) ToMeter() float64 {
	return p.to_meter
}
func (p *pj) FromGreenwich() float64 {
	return p.from_greenwich
}
func (p *pj) Radius() float64 {
	return p.a
}
