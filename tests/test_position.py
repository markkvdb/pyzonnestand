"""Check that the sun position calculations are correct."""

from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

from pyzonnestand.position import (
    _abberation_correction,
    _calendar_time,
    _geocentric_position,
    _greenwich_sidereal_time,
    _heliocentric_position,
    _julian_century,
    _julian_day,
    _julian_ephemeris_day,
    _julian_millennium,
    _nutation_obliquity,
    _sun_longitude,
    _sun_ra_decl,
    _sun_topo_azimuth_zenith,
    _sun_topo_ra_decl_hour,
    sunpos,
)


def test_position() -> None:
    test_file = Path(__file__).parent / "data" / "sample.txt"
    # Parse and compare results from https://midcdmz.nrel.gov/solpos/spa.html
    param_names = [
        "syear",
        "smonth",
        "sday",
        "eyear",
        "emonth",
        "eday",
        "otype",
        "step",
        "stepunit",
        "hr",
        "min",
        "sec",
        "latitude",
        "longitude",
        "timezone",
        "elev",
        "press",
        "temp",
        "dut1",
        "deltat",
        "azmrot",
        "slope",
        "refract",
    ]
    param_dtype = np.dtype([(name, float) for name in param_names])
    params = np.loadtxt(test_file, param_dtype, delimiter=",", skiprows=2, max_rows=1)

    row_type = np.dtype(
        [
            ("Date_M/D/YYYY", "S10"),
            ("Time_H:MM:SS", "S8"),
            ("Topo_zen", float),
            ("Topo_az", float),
            ("Julian_day", float),
            ("Julian_century", float),
            ("Julian_ephemeris_day", float),
            ("Julian_ephemeris_century", float),
            ("Julian_ephemeris_millennium", float),
            ("Earth_heliocentric_longitude", float),
            ("Earth_heliocentric_latitude", float),
            ("Earth_radius_vector", float),
            ("Geocentric_longitude", float),
            ("Geocentric_latitude", float),
            ("Mean_elongation", float),
            ("Mean_anomaly_sun", float),
            ("Mean_anomaly_moon", float),
            ("Argument_latitude_moon", float),
            ("Ascending_longitude_moon", float),
            ("Nutation_longitude", float),
            ("Nutation_obliquity", float),
            ("Ecliptic_mean_obliquity", float),
            ("Ecliptic_true_obliquity", float),
            ("Aberration_correction", float),
            ("Apparent_sun_longitude", float),
            ("Greenwich_mean_sidereal_time", float),
            ("Greenwich_sidereal_time", float),
            ("Geocentric_sun_right_ascension", float),
            ("Geocentric_sun_declination", float),
            ("Observer_hour_angle", float),
            ("Sun_equatorial_horizontal_parallax", float),
            ("Sun_right_ascension_parallax", float),
            ("Topo_sun_declination", float),
            ("Topo_sun_right_ascension", float),
            ("Topo_local_hour_angle", float),
            ("Topo_elevation_angle_uncorrected", float),
            ("Atmospheric_refraction_correction", float),
            ("Topo_elevation_angle_corrected", float),
            ("Equation_of_time", float),
            ("Sunrise_hour_angle", float),
            ("Sunset_hour_angle", float),
            ("Sun_transit_altitude", float),
        ]
    )

    true_data = np.loadtxt(test_file, row_type, delimiter=",", skiprows=4)

    def to_datetime(date_time_pair):
        s = str(b" ".join(date_time_pair), "UTF-8")
        return datetime.strptime(s, "%m/%d/%Y %H:%M:%S")

    def angle_diff(a1, a2, period=2 * np.pi):
        """(a1 - a2 + d) % (2*d) - d; d = period/2"""
        d = period / 2
        return ((a1 - a2 + d) % (period)) - d

    dts = [
        to_datetime(dt_pair) for dt_pair in true_data[["Date_M/D/YYYY", "Time_H:MM:SS"]]
    ]
    lat, lon, elev, temp, press, deltat = (
        params["latitude"],
        params["longitude"],
        params["elev"],
        params["temp"],
        params["press"],
        params["deltat"],
    )
    all_errs = []
    for dt, truth in zip(dts, true_data):
        t = _calendar_time(dt)
        jd = _julian_day(t)  # Julian_day
        jde = _julian_ephemeris_day(jd, deltat)  # Julian_ephemeris_day
        jce = _julian_century(jde)  # Julian_ephemeris_century
        jme = _julian_millennium(jce)  # Julian_ephemeris_millenium
        L, B, R = _heliocentric_position(
            jme
        )  # Earth_heliocentric_longitude, Earth_heliocentric_latitude, Earth_radius_vector
        delta_psi, epsilon = _nutation_obliquity(
            jce
        )  # Nutation_longitude, Ecliptic_true_obliquity
        theta, beta = _geocentric_position(
            (L, B, R)
        )  # Geocentric_longitude, Geocentric_latitude
        delta_tau = _abberation_correction(R)  # Aberration_correction
        llambda, beta = _sun_longitude(
            (L, B, R), delta_psi
        )  # Apparent_sun_longitude, Geocentric_latitude (identical to previous)
        v = _greenwich_sidereal_time(jd, delta_psi, epsilon)  # Greenwich_sidereal_time
        alpha, delta = _sun_ra_decl(
            llambda, epsilon, beta
        )  # Geocentric_sun_right_ascension, Geocentric_sun_declination
        alpha_p, delta_p, H_p = _sun_topo_ra_decl_hour(
            lat, lon, elev, jd, deltat
        )  # Topo_sun_right_ascension, Topo_sun_declination, Topo_local_hour_angle
        az, zen, delta_e = _sun_topo_azimuth_zenith(
            lat, delta_p, H_p, temp, press
        )  # Topo_az, Topo_zen, Atmospheric_refraction_correction

        jd_err = jd - truth["Julian_day"]
        jde_err = jde - truth["Julian_ephemeris_day"]
        jce_err = jce - truth["Julian_ephemeris_century"]
        jme_err = jme - truth["Julian_ephemeris_millennium"]
        L_err = L - truth["Earth_heliocentric_longitude"]
        B_err = B - truth["Earth_heliocentric_latitude"]
        R_err = R - truth["Earth_radius_vector"]
        delta_psi_err = delta_psi - truth["Nutation_longitude"]
        epsilon_err = epsilon - truth["Ecliptic_true_obliquity"]
        theta_err = theta - truth["Geocentric_longitude"]
        beta_err = beta - truth["Geocentric_latitude"]
        delta_tau_err = delta_tau - truth["Aberration_correction"]
        lambda_err = llambda - truth["Apparent_sun_longitude"]
        v_err = v - truth["Greenwich_sidereal_time"]
        alpha_err = alpha - truth["Geocentric_sun_right_ascension"]
        delta_err = delta - truth["Geocentric_sun_declination"]
        alpha_prime_err = alpha_p - truth["Topo_sun_right_ascension"]
        delta_prime_err = delta_p - truth["Topo_sun_declination"]
        H_prime_err = angle_diff(H_p, truth["Topo_local_hour_angle"], 360)
        az_err = angle_diff(az, truth["Topo_az"], 360)
        delta_e_err = delta_e - truth["Atmospheric_refraction_correction"]
        zen_err = zen - truth["Topo_zen"]
        all_errs.append(
            [
                jd_err,
                jde_err,
                jce_err,
                jme_err,
                L_err,
                B_err,
                R_err,
                delta_psi_err,
                epsilon_err,
                theta_err,
                beta_err,
                delta_tau_err,
                lambda_err,
                v_err,
                alpha_err,
                delta_err,
                alpha_prime_err,
                delta_prime_err,
                H_prime_err,
                az_err,
                delta_e_err,
                zen_err,
            ]
        )
    rms_err = np.sqrt(np.mean(np.array(all_errs) ** 2, 0))
    err_names = [
        "Julian day",
        "Julian ephemeris day",
        "Julian ephemeris century",
        "Julian ephemeris millennium",
        "Earth heliocentric longitude",
        "Earth heliocentric latitude",
        "Earth radius vector",
        "Nutation longitude",
        "Ecliptic true obliquity",
        "Geocentric longitude",
        "Geocentric latitude",
        "Aberration correction",
        "Apparent sun longitude",
        "Greenwich sidereal time",
        "Geocentric sun right ascension",
        "Geocentric sun declination",
        "Topo sun right ascension",
        "Topo sun declination",
        "Topo local hour angle",
        "Topo az",
        "Atmospheric_refraction_correction",
        "Topo zen",
    ]
    print("RMS Errors")
    for n, e in zip(err_names, rms_err):
        assert e < 1e-4, f"{n}: {e}"


def test_sunpos() -> None:
    """Test sun position calculation."""
    start = datetime(2019, 1, 1, 12, 0, 0)
    dt = [start + timedelta(hours=i) for i in range(48)]
    x = sunpos(
        dt=dt,
        latitude=52.0,
        longitude=5.0,
        elevation=0.0,
    )

    assert np.array(x).shape == (5, 48)
