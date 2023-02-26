#!/usr/bin/env python3
import json
import re
from datetime import datetime, timedelta
from enum import Enum
from math import radians
from typing import Any, Dict, Iterator, List, NamedTuple, Optional, Tuple

astronomical_unit = 149597870700.
century = 36525 * 86400

dwarf_planets = (
    'Ceres',
    'Orcus',
    'Pluto',
    'Haumea',
    'Quaoar',
    'Makemake',
    'Gonggong',
    'Eris',
    'Sedna',
)


class PlanetPhysParDataRow(NamedTuple):
    name: str
    equatorial_radius: float
    mean_radius: float
    mass: float
    bulk_density: float
    sidereal_rotation_period: float
    sidereal_orbital_period: float
    magnitude: float
    geometric_albedo: float
    equatorial_gravity: float
    escape_velocity: float


class FactsheetDataRow(NamedTuple):
    name: str
    # TODO: many other values are available
    # TODO: precession
    # north pole of rotation
    north_pole_right_ascension: Optional[float]
    north_pole_declination: Optional[float]


class PElemDataRow(NamedTuple):
    name: str
    semi_major_axis: Tuple[float, float]
    eccentricity: Tuple[float, float]
    inclination: Tuple[float, float]
    mean_longitude: Tuple[float, float]
    longitude_of_periapsis: Tuple[float, float]
    longitude_of_ascending_node: Tuple[float, float]


class SatPhysParDataRow(NamedTuple):
    name: str
    primary: str
    gravitational_parameter: float
    mean_radius: float
    mean_density: Optional[float]


class ReferencePlane(str, Enum):
    ECLIPTIC = 'ECLIPTIC'
    EQUATORIAL = 'EQUATORIAL'
    LAPLACE = 'LAPLACE'


class SatElemDataRow(NamedTuple):
    name: str
    primary: str
    reference_plane: ReferencePlane
    epoch: datetime
    semi_major_axis: float
    eccentricity: float
    argument_of_periapsis: float
    mean_anomaly_at_epoch: float
    inclination: float
    longitude_of_ascending_node: float
    period: float
    argument_of_periapsis_precession_period: float
    longitude_of_ascending_node_precession_period: float
    north_pole_right_ascension: Optional[float]
    north_pole_declination: Optional[float]
    tilt: Optional[float]


class SBDBDataRow(NamedTuple):
    name: str
    primary: str
    epoch: datetime
    eccentricity: float
    semi_major_axis: float
    periapsis: float
    inclination: float
    longitude_of_ascending_node: float
    argument_of_periapsis: float
    mean_anomaly_at_epoch: float
    orbital_period: float
    mean_motion: float
    apoapsis: float

    absolute_magnitude: Optional[float] = None
    magnitude_slope: Optional[float] = None
    radius: Optional[float] = None
    gravitational_parameter: Optional[float] = None
    bulk_density: Optional[float] = None
    rotation_period: Optional[float] = None
    north_pole_right_ascension: Optional[float] = None
    north_pole_declination: Optional[float] = None
    geometric_albedo: Optional[float] = None


class SheppardMoonElements(NamedTuple):
    number: str
    name: Optional[str]
    designation: Optional[str]
    semi_major_axis: float
    inclination: float
    eccentricity: float
    argument_of_periapsis: float
    longitude_of_ascending_node: float
    mean_anomaly: float
    period: float
    magnitude: float
    radius: float
    year: int


def planet_phys_par_G() -> float:
    with open('data/planet_phys_par.html') as f:
        html = f.read()

    match = re.search(r'(?s).*G=(6\.674[0-9]*)', html)
    assert match is not None
    return float(match.group(1)) * 1e-11


def parse_planet_phys_par() -> Iterator[PlanetPhysParDataRow]:
    with open('data/planet_phys_par.html') as f:
        html = f.read()

    # iterate over rows

    matches = re.findall(r'(?s)<tr>\s*<td><b>([A-Z]\w*)</b></td>(.*?)</tr>', html)
    assert len(matches) >= 8
    for planet, data in matches:
        # extract data cells from row (ingore precision information)
        matches = re.findall(r'(?s)<td[^>]*>\s*(\S+)', data)

        # unpack values
        (
            equatorial_radius, mean_radius, mass, bulk_density,
            sidereal_rotation_period, sidereal_orbital_period, magnitude, albedo,
            equatorial_gravity, escape_velocity,
        ) = [float(match) for match in matches]

        if planet == 'Pluto':
            # see IAU definition of positive pole for dwarf planets
            sidereal_rotation_period = abs(sidereal_rotation_period)

        # normalize to SI units and export
        yield PlanetPhysParDataRow(
            name=planet,
            equatorial_radius=equatorial_radius * 1e3,
            mean_radius=mean_radius * 1e3,
            mass=mass * 1e24,
            bulk_density=bulk_density * 1e3,
            sidereal_rotation_period=sidereal_rotation_period * 86400,
            sidereal_orbital_period=sidereal_orbital_period * 86400 * 365.25,
            magnitude=magnitude,
            geometric_albedo=albedo,
            equatorial_gravity=equatorial_radius,
            escape_velocity=escape_velocity * 1e3,
        )


def parse_factsheet(filename: str) -> FactsheetDataRow:
    with open(f'data/factsheet/{filename}') as f:
        html = f.read()

    match = re.search(r'<title>([^<]*) Fact Sheet</title>', html)
    assert match is not None
    name = match.group(1)

    # orientation of north pole (axial tilt)

    # extract right ascension
    match = re.search(r'Right Ascension *: *([\-0-9\.]+)', html)
    if match is None:
        north_pole_right_ascension = None
    else:
        north_pole_right_ascension = radians(float(match.group(1)))

    # extract declination
    match = re.search(r'Declination *: *([\-0-9\.]+)', html)
    if match is None:
        north_pole_declination = None
    else:
        north_pole_declination = radians(float(match.group(1)))

    return FactsheetDataRow(
        name=name,
        north_pole_right_ascension=north_pole_right_ascension,
        north_pole_declination=north_pole_declination,
    )


def parse_factsheets() -> Iterator[FactsheetDataRow]:
    with open('data/factsheet/index.html') as f:
        html = f.read()

    matches = re.findall(r'href="([^"/]*fact.html)"', html)
    assert matches
    for filename in set(matches):
        yield parse_factsheet(filename)


def parse_p_elem(lines: List[str]) -> Iterator[PElemDataRow]:
    # iterate lines that represent rows in the table
    for i in range(0, len(lines), 2):
        # each row is made of two lines
        first_line = lines[i]
        second_line = lines[i + 1]

        # names can contains space ("EM Bary")
        name, *str_elements = re.split(r'\s{2,}', first_line)
        str_changes = second_line.split()

        # parse values
        elements = [float(element) for element in str_elements]
        changes = [float(change) for change in str_changes]

        # convert denominator of rates of change from century to second
        changes = [change / century for change in changes]

        # unpack
        (
            semi_major_axis, eccentricity, inclination, mean_longitude,
            longitude_of_periapsis, longitude_of_ascending_node,
        ) = elements
        (
            r_semi_major_axis, r_eccentricity, r_inclination, r_mean_longitude,
            r_longitude_of_periapsis, r_longitude_of_ascending_node,
        ) = changes

        # normalize to SI units
        semi_major_axis *= astronomical_unit
        r_semi_major_axis *= astronomical_unit

        inclination = radians(inclination)
        r_inclination = radians(r_inclination)

        mean_longitude = radians(mean_longitude)
        r_mean_longitude = radians(r_mean_longitude)

        longitude_of_periapsis = radians(longitude_of_periapsis)
        r_longitude_of_periapsis = radians(r_longitude_of_periapsis)

        longitude_of_ascending_node = radians(longitude_of_ascending_node)
        r_longitude_of_ascending_node = radians(r_longitude_of_ascending_node)

        # export
        yield PElemDataRow(
            name=name,
            semi_major_axis=(semi_major_axis, r_semi_major_axis),
            eccentricity=(eccentricity, r_eccentricity),
            inclination=(inclination, r_inclination),
            mean_longitude=(mean_longitude, mean_longitude),
            longitude_of_periapsis=(longitude_of_periapsis, r_longitude_of_periapsis),
            longitude_of_ascending_node=(longitude_of_ascending_node, r_longitude_of_ascending_node),
        )


def parse_p_elem_t1() -> Iterator[PElemDataRow]:
    with open('data/approx_pos.html') as f:
        lines = f.readlines()
    start = lines.index('<pre> \n') + 4
    stop = lines.index('</pre>\n') - 2
    yield from parse_p_elem(lines[start:stop])


def parse_sat_phys_par() -> Iterator[SatPhysParDataRow]:
    with open('data/sat_phys_par.html') as f:
        html = f.read()

    matches = re.findall(r'(?s)<tr>\s*<td class="text-left">([A-Z]\w*)</td>(.*?)</tr>', html)
    assert len(matches) > 10
    for primary, data in matches:
        values = re.findall(r'(?s)<td[^>]*>\s*([^<]+?)\s*<', data)
        name, code, gravitational_parameter, GM_ref, mean_radius, _, mean_density, _ = values

        # normalize to SI units and export
        yield SatPhysParDataRow(
            name=name,
            primary=primary,
            gravitational_parameter=float(gravitational_parameter) * 1e9,
            mean_radius=float(mean_radius) * 1e3,
            mean_density=None if mean_density == 'n/a' else float(mean_density) * 1e3,
        )


def parse_sat_elem() -> Iterator[SatElemDataRow]:
    with open('data/sat_elem.html') as f:
        html = f.read()

    matches = re.findall(r'(?s)<tr>\s*<td>\d+</td>(.*?)</tr>', html)
    assert len(matches) > 10
    for data in matches:
        values = re.findall(r'(?s)<td[^>]*>\s*([^<]*?)\s*<', data)
        primary, name, code, _, frame, epoch, *elements, ref = values

        if 'ecliptic' in data:
            reference_plane = ReferencePlane.ECLIPTIC
        elif 'Laplace' in data:
            reference_plane = ReferencePlane.LAPLACE
        elif 'equatorial' in data:
            # NOTE: the elements of Pluto's moons are barycentric except
            # for Charon, for which they are plateocentric
            reference_plane = ReferencePlane.EQUATORIAL
        else:
            assert False

        # parse epoch
        date, fraction = epoch.split('.')
        epoch = datetime.strptime(date, '%Y-%m-%d') + timedelta(days=int(fraction) / 10)

        # unpack
        (
            semi_major_axis,
            eccentricity,
            argument_of_periapsis,
            mean_anomaly_at_epoch,
            inclination,
            longitude_of_ascending_node,
            period,
            argument_of_periapsis_precession_period,
            longitude_of_ascending_node_precession_period,
            *laplace_elements,
        ) = [0 if element in ('', '-') else float(element) for element in elements]

        # special case for Laplace reference plane
        north_pole_right_ascension: Optional[float]
        north_pole_declination: Optional[float]
        tilt: Optional[float]
        if reference_plane == 'laplace':
            laplace_elements = [radians(element) for element in laplace_elements]
            north_pole_right_ascension, north_pole_declination, tilt = laplace_elements
        else:
            north_pole_right_ascension = None
            north_pole_declination = None
            tilt = None

        # normalize to SI units and export
        yield SatElemDataRow(
            name=name,
            primary=primary,
            reference_plane=reference_plane,
            epoch=epoch,
            semi_major_axis=semi_major_axis * 1e3,
            eccentricity=eccentricity,
            argument_of_periapsis=radians(argument_of_periapsis),
            mean_anomaly_at_epoch=radians(mean_anomaly_at_epoch),
            inclination=radians(inclination),
            longitude_of_ascending_node=radians(longitude_of_ascending_node),
            period=period * 86400,
            argument_of_periapsis_precession_period=(
                argument_of_periapsis_precession_period * 365.25 * 86400
            ),
            longitude_of_ascending_node_precession_period=(
                longitude_of_ascending_node_precession_period * 365.25 * 86400
            ),
            north_pole_right_ascension=north_pole_right_ascension,
            north_pole_declination=north_pole_declination,
            tilt=tilt,
        )


def parse_sbdb(name: str) -> SBDBDataRow:
    with open(f'data/sbdb/{name}.json') as f:
        data = json.load(f)

    parameters: Dict = {
        'name': name,
        'primary': 'Sun',
    }

    # epoch
    julian_epoch = float(data['orbit']['epoch'])
    julian_J2000 = 2451545.0
    J2000 = datetime(2000, 1, 1, 12)
    parameters['epoch'] = J2000 + timedelta(days=julian_epoch - julian_J2000)

    # orbital elements
    for element in data['orbit']['elements']:
        pname = element['name']
        value = element['value']
        if pname == 'e':
            parameters['eccentricity'] = float(value)
        elif pname == 'a':
            parameters['semi_major_axis'] = float(value) * astronomical_unit
        elif pname == 'q':
            parameters['periapsis'] = float(value) * astronomical_unit
        elif pname == 'i':
            parameters['inclination'] = radians(float(value))
        elif pname == 'om':
            parameters['longitude_of_ascending_node'] = radians(float(value))
        elif pname == 'w':
            parameters['argument_of_periapsis'] = radians(float(value))
        elif pname == 'ma':
            parameters['mean_anomaly_at_epoch'] = radians(float(value))
        elif pname == 'tp':  # t_p, time to periapsis
            pass
        elif pname == 'per':
            parameters['orbital_period'] = float(value)
        elif pname == 'n':
            parameters['mean_motion'] = radians(float(value)) / 86400
        elif pname == 'ad':
            parameters['apoapsis'] = float(value) * astronomical_unit
        else:
            print(f'Parsing SBDB orbital elements for {name}, dit not expect "{pname}" ({value})')

    for parameter in data['phys_par']:
        pname = parameter['name']
        value = parameter['value']
        if pname == 'H':
            parameters['absolute_magnitude'] = float(value)
        elif pname == 'G':
            parameters['magnitude_slope'] = float(value)
        elif pname == 'diameter':
            parameters['radius'] = float(value) * 1e3 / 2.
        elif pname == 'extent':
            pass
        elif pname == 'GM':
            parameters['gravitational_parameter'] = float(value) * 1e9
        elif pname == 'density':
            parameters['bulk_density'] = float(value) * 1e3
        elif pname == 'rot_per':
            parameters['rotation_period'] = float(value) * 3600.
        elif pname == 'pole':
            ra, dec = value.split('/')
            parameters['north_pole_right_ascension'] = float(ra)
            parameters['north_pole_declination'] = float(dec)
        elif pname == 'albedo':
            parameters['geometric_albedo'] = float(value)
        elif pname == 'BV':
            pass
        elif pname == 'UB':
            pass
        elif pname == 'spec_T':
            pass
        elif pname == 'spec_B':
            pass
        else:
            print(f'Parsing SBDB physical parameters for {name}, dit not expect "{pname}" ({value})')

    return SBDBDataRow(**parameters)


def parse_sheppard_moons(primary: str) -> Iterator[SheppardMoonElements]:
    with open(f'data/{primary.lower()}moons') as f:
        html = f.read()

    # fields: Number, Name, Designation, a, i, e, Peri, Node, M, Period, mag, Size, year
    # NOTE: this assumes there is no space in the name
    # NOTE: “S” in “S/2000 J1” stands for “Satellite”
    pattern = r'&gt;(?:([IVXLCDM\?]+|P\d) +)?(?:(\S*) +)?(?:(S/\d{4} .\d+) +)?(\d+)' + ' +([0-9.]+)' * 9
    matches = re.findall(pattern, html)
    assert matches
    for number, name, designation, a, i, e, om, Om, M, T, mag, size, year in matches:
        # put back space in designation ("S/2000 J01" → "S/2000 J 1")
        if designation:
            left, right = designation.split()
            designation = f'{left} {right[0]} {right[1:].lstrip("0")}'
        yield SheppardMoonElements(
            number=number,
            name=name,
            designation=designation,
            semi_major_axis=float(a) * 1e3,
            inclination=radians(float(i)),
            eccentricity=float(e),
            argument_of_periapsis=radians(float(om)),
            longitude_of_ascending_node=radians(float(Om)),
            mean_anomaly=radians(float(M)),
            period=float(T) * 86400,
            magnitude=float(mag),
            radius=float(size) * 500,
            year=int(year),
        )


def main() -> None:
    row: Any = None
    print(planet_phys_par_G())

    for row in parse_planet_phys_par():
        print(row)

    for row in parse_factsheets():
        print(row)

    for row in parse_p_elem_t1():
        print(row)

    for row in parse_sat_phys_par():
        print(row)

    for row in parse_sat_elem():
        print(row)

    for planet in ('Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'):
        for moon in parse_sheppard_moons(planet):
            print(moon)

    for dwarf_planet in dwarf_planets:
        print(parse_sbdb(dwarf_planet))


if __name__ == '__main__':
    main()
