// Copyright 2015 Sam L'ecuyer. All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

package projectron

import (
	"strconv"
	"strings"
	"errors"
	"math"
)

type translator func(float64, float64) (float64, float64, error)

const (
	sixth   float64 = 1 / 6
	ra4             = 17 / 360
	ra6             = 67 / 3024
	rv4             = 5 / 72
	rv6             = 55 / 1296
	sec2rad         = 4.84813681109535993589914102357e-6
	epsln           = 1.0e-10
)

type datumType int

const (
	PJD_UNKNOWN datumType = 0
	PJD_GRIDSHIFT
	PJD_7PARAM
	PJD_3PARAM
	PJD_WGS84
)

var ErrUnsupportedProj = errors.New("This is not a supported Projection")
var ErrUnknownDatum = errors.New("This is not a supported datum")
var ErrInvalidParam = errors.New("We encountered an illegal parameter")

var hugeVal = math.Inf(1)

type paramset map[string]string

func (p paramset) string(s string) (v string, ok bool) {
	v, ok = p[s]
	return v, ok
}
func (p paramset) bool(s string) (b bool, okay bool) {
	var err error
	if v, ok := p[s]; ok {
		if v == "" {
			return true, true
		}
		b, err = strconv.ParseBool(v)
		okay = err == nil
	}
	return
}
func (p paramset) float(s string) (f float64, okay bool) {
	var err error
	if v, ok := p[s]; ok {
		f, err = strconv.ParseFloat(v, 64)
		okay = err == nil
	}
	return
}
func (p paramset) degree(s string) (f float64, okay bool) {
	// var err error
	if v, ok := p[s]; ok {
		return parseDegreeString(v) * d2r, ok
	}
	return
}

func parseDegreeString(ds string) float64 {
	var res float64
	idx := strings.Index(ds, "d")
	if idx >= 0 {
		f, _ := strconv.ParseFloat(ds[0:idx], 64)
		res += f
		ds = ds[idx+1:]
	} else {
		res, _ = strconv.ParseFloat(ds, 64)
	}
	idx = strings.Index(ds, "'")
	if idx >= 0 {
		f, _ := strconv.ParseFloat(ds[0:idx], 64)
		res += f / 60
		ds = ds[idx+1:]
	}
	idx = strings.Index(ds, "\"")
	if idx >= 0 {
		f, _ := strconv.ParseFloat(ds[0:idx], 64)
		res += f / 360
		ds = ds[idx+1:]
	}
	if strings.HasSuffix(ds, "W") || strings.HasSuffix(ds, "S") {
		res *= -1
	}
	return res
}

func keyVal(s string) (key string, val string) {
	defs := strings.Split(s, "=")
	key = defs[0]
	if len(defs) == 2 {
		val = defs[1]
	}
	return
}

type datum struct {
	id, definition, ellipse, comments string
}

// https://github.com/OSGeo/proj.4/blob/6eabf3c854f3786fbaa1b645ea9e9aae6e4522e4/src/pj_datums.c#L41-L62
var datums_list = map[string]datum{
	"WGS84": datum{"WGS84", "towgs84=0,0,0", "WGS84", ""},
	"GGRS87": datum{"GGRS87", "towgs84=-199.87,74.79,246.62", "GRS80",
		"Greek_Geodetic_Reference_System_1987"},
	"NAD83": datum{"NAD83", "towgs84=0,0,0", "GRS80",
		"North_American_Datum_1983"},
	"NAD27": datum{"NAD27", "nadgrids=@conus,@alaska,@ntv2_0.gsb,@ntv1_can.dat",
		"clrk66",
		"North_American_Datum_1927"},
	"potsdam": datum{"potsdam", "towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7",
		"bessel",
		"Potsdam Rauenberg 1950 DHDN"},
	"carthage": datum{"carthage", "towgs84=-263.0,6.0,431.0", "clrk80ign",
		"Carthage 1934 Tunisia"},
	"hermannskogel": datum{"hermannskogel", "towgs84=577.326,90.129,463.919,5.137,1.474,5.297,2.4232",
		"bessel",
		"Hermannskogel"},
	"ire65": datum{"ire65", "towgs84=482.530,-130.596,564.557,-1.042,-0.214,-0.631,8.15",
		"mod_airy",
		"Ireland 1965"},
	"nzgd49": datum{"nzgd49", "towgs84=59.47,-5.04,187.44,0.47,-0.1,1.024,-4.5993",
		"intl", "New Zealand Geodetic Datum 1949"},
	"OSGB36": datum{"OSGB36", "towgs84=446.448,-125.157,542.060,0.1502,0.2470,0.8421,-20.4894",
		"airy", "Airy 1830"},
}

type ellipse struct {
	id, major, ell, name string
}

// https://github.com/OSGeo/proj.4/blob/6eabf3c854f3786fbaa1b645ea9e9aae6e4522e4/src/pj_ellps.c#L7-L49
var ellipse_list = map[string]ellipse{
	"MERIT":     ellipse{"MERIT", "a=6378137.0", "rf=298.257", "MERIT 1983"},
	"SGS85":     ellipse{"SGS85", "a=6378136.0", "rf=298.257", "Soviet Geodetic System 85"},
	"GRS80":     ellipse{"GRS80", "a=6378137.0", "rf=298.257222101", "GRS 1980(IUGG, 1980)"},
	"IAU76":     ellipse{"IAU76", "a=6378140.0", "rf=298.257", "IAU 1976"},
	"airy":      ellipse{"airy", "a=6377563.396", "b=6356256.910", "Airy 1830"},
	"APL4.9":    ellipse{"APL4.9", "a=6378137.0.", "rf=298.25", "Appl. Physics. 1965"},
	"NWL9D":     ellipse{"NWL9D", "a=6378145.0.", "rf=298.25", "Naval Weapons Lab., 1965"},
	"mod_airy":  ellipse{"mod_airy", "a=6377340.189", "b=6356034.446", "Modified Airy"},
	"andrae":    ellipse{"andrae", "a=6377104.43", "rf=300.0", "Andrae 1876 (Den., Iclnd.)"},
	"aust_SA":   ellipse{"aust_SA", "a=6378160.0", "rf=298.25", "Australian Natl & S. Amer. 1969"},
	"GRS67":     ellipse{"GRS67", "a=6378160.0", "rf=298.2471674270", "GRS 67(IUGG 1967)"},
	"bessel":    ellipse{"bessel", "a=6377397.155", "rf=299.1528128", "Bessel 1841"},
	"bess_nam":  ellipse{"bess_nam", "a=6377483.865", "rf=299.1528128", "Bessel 1841 (Namibia)"},
	"clrk66":    ellipse{"clrk66", "a=6378206.4", "b=6356583.8", "Clarke 1866"},
	"clrk80":    ellipse{"clrk80", "a=6378249.145", "rf=293.4663", "Clarke 1880 mod."},
	"clrk80ign": ellipse{"clrk80ign", "a=6378249.2", "rf=293.4660212936269", "Clarke 1880 (IGN)."},
	"CPM":       ellipse{"CPM", "a=6375738.7", "rf=334.29", "Comm. des Poids et Mesures 1799"},
	"delmbr":    ellipse{"delmbr", "a=6376428.", "rf=311.5", "Delambre 1810 (Belgium)"},
	"engelis":   ellipse{"engelis", "a=6378136.05", "rf=298.2566", "Engelis 1985"},
	"evrst30":   ellipse{"evrst30", "a=6377276.345", "rf=300.8017", "Everest 1830"},
	"evrst48":   ellipse{"evrst48", "a=6377304.063", "rf=300.8017", "Everest 1948"},
	"evrst56":   ellipse{"evrst56", "a=6377301.243", "rf=300.8017", "Everest 1956"},
	"evrst69":   ellipse{"evrst69", "a=6377295.664", "rf=300.8017", "Everest 1969"},
	"evrstSS":   ellipse{"evrstSS", "a=6377298.556", "rf=300.8017", "Everest (Sabah & Sarawak)"},
	"fschr60":   ellipse{"fschr60", "a=6378166.", "rf=298.3", "Fischer (Mercury Datum) 1960"},
	"fschr60m":  ellipse{"fschr60m", "a=6378155.", "rf=298.3", "Modified Fischer 1960"},
	"fschr68":   ellipse{"fschr68", "a=6378150.", "rf=298.3", "Fischer 1968"},
	"helmert":   ellipse{"helmert", "a=6378200.", "rf=298.3", "Helmert 1906"},
	"hough":     ellipse{"hough", "a=6378270.0", "rf=297.", "Hough"},
	"intl":      ellipse{"intl", "a=6378388.0", "rf=297.", "International 1909 (Hayford)"},
	"krass":     ellipse{"krass", "a=6378245.0", "rf=298.3", "Krassovsky, 1942"},
	"kaula":     ellipse{"kaula", "a=6378163.", "rf=298.24", "Kaula 1961"},
	"lerch":     ellipse{"lerch", "a=6378139.", "rf=298.257", "Lerch 1979"},
	"mprts":     ellipse{"mprts", "a=6397300.", "rf=191.", "Maupertius 1738"},
	"new_intl":  ellipse{"new_intl", "a=6378157.5", "b=6356772.2", "New International 1967"},
	"plessis":   ellipse{"plessis", "a=6376523.", "b=6355863.", "Plessis 1817 (France)"},
	"SEasia":    ellipse{"SEasia", "a=6378155.0", "b=6356773.3205", "Southeast Asia"},
	"walbeck":   ellipse{"walbeck", "a=6376896.0", "b=6355834.8467", "Walbeck"},
	"WGS60":     ellipse{"WGS60", "a=6378165.0", "rf=298.3", "WGS 60"},
	"WGS66":     ellipse{"WGS66", "a=6378145.0", "rf=298.25", "WGS 66"},
	"WGS72":     ellipse{"WGS72", "a=6378135.0", "rf=298.26", "WGS 72"},
	"WGS84":     ellipse{"WGS84", "a=6378137.0", "rf=298.257223563", "WGS 84"},
	"sphere":    ellipse{"sphere", "a=6370997.0", "b=6370997.0", "Normal Sphere (r=6370997)"},
}

type unit struct {
	id       string
	to_meter float64
	name     string
}

// https://github.com/OSGeo/proj.4/blob/6eabf3c854f3786fbaa1b645ea9e9aae6e4522e4/src/pj_units.c#L9-L30
var units_list = map[string]unit{
	"km":     unit{"km", 1000, "Kilometer"},
	"m":      unit{"m", 1.0, "Meter"},
	"dm":     unit{"dm", 0.1, "Decimeter"},
	"cm":     unit{"cm", 0.01, "Centimeter"},
	"mm":     unit{"mm", 0.001, "Millimeter"},
	"kmi":    unit{"kmi", 1852.0, "International Nautical Mile"},
	"in":     unit{"in", 0.0254, "International Inch"},
	"ft":     unit{"ft", 0.3048, "International Foot"},
	"yd":     unit{"yd", 0.9144, "International Yard"},
	"mi":     unit{"mi", 1609.344, "International Statute Mile"},
	"fath":   unit{"fath", 1.8288, "International Fathom"},
	"ch":     unit{"ch", 20.1168, "International Chain"},
	"link":   unit{"link", 0.201168, "International Link"},
	"us-in":  unit{"us-in", 0.0254000508, "U.S. Surveyor's Inch"},
	"us-ft":  unit{"us-ft", 0.304800609601219, "U.S. Surveyor's Foot"},
	"us-yd":  unit{"us-yd", 0.914401828803658, "U.S. Surveyor's Yard"},
	"us-ch":  unit{"us-ch", 20.11684023368047, "U.S. Surveyor's Chain"},
	"us-mi":  unit{"us-mi", 1609.347218694437, "U.S. Surveyor's Statute Mile"},
	"ind-yd": unit{"ind-yd", 0.91439523, "Indian Yard"},
	"ind-ft": unit{"ind-ft", 0.30479841, "Indian Foot"},
	"ind-ch": unit{"ind-ch", 20.11669506, "Indian Chain"},
}

type prime_meridian struct {
	id, defn string
}

// https://github.com/OSGeo/proj.4/blob/6eabf3c854f3786fbaa1b645ea9e9aae6e4522e4/src/pj_units.c#L9-L30
var pm_list = map[string]prime_meridian{
	"greenwich": prime_meridian{"greenwich", "0dE"},
	"lisbon":    prime_meridian{"lisbon", "9d07'54.862\"W"},
	"paris":     prime_meridian{"paris", "2d20'14.025\"E"},
	"bogota":    prime_meridian{"bogota", "74d04'51.3\"W"},
	"madrid":    prime_meridian{"madrid", "3d41'16.58\"W"},
	"rome":      prime_meridian{"rome", "12d27'8.4\"E"},
	"bern":      prime_meridian{"bern", "7d26'22.5\"E"},
	"jakarta":   prime_meridian{"jakarta", "106d48'27.79\"E"},
	"ferro":     prime_meridian{"ferro", "17d40'W"},
	"brussels":  prime_meridian{"brussels", "4d22'4.71\"E"},
	"stockholm": prime_meridian{"stockholm", "18d3'29.8\"E"},
	"athens":    prime_meridian{"athens", "23d42'58.815\"E"},
	"oslo":      prime_meridian{"oslo", "10d43'22.5\"E"},
}
