"""
Nancy Grace Roman Space Telescope - Ground Telemetry Workstation (v3 Beta 2)
Lightweight Python telemetry viewer, Lupton RGB composite engine, Dual Pixel
tracker, and Plain-English Science Translator for MAST WFI exposures.
"""

import os
import gc
import math
import smtplib
import threading
import urllib.request
import json
from datetime import datetime, timezone
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.utils import formatdate, make_msgid

import tkinter as tk
from PIL import Image, ImageTk

from astropy.io import fits
from astropy.coordinates import EarthLocation, SkyCoord, AltAz
from astropy.time import Time
from astropy.visualization import ZScaleInterval, AsinhStretch, ImageNormalize, make_lupton_rgb
import astropy.units as u
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# =====================================================================
# STATION TELEMETRY & EMAIL CONFIGURATION (Set via Environment Variables)
# =====================================================================
STATION_NAME = os.getenv("ROMAN_STATION_NAME", "ROMAN GROUND COMMAND // GROUND STATION")
STATION_LAT = float(os.getenv("ROMAN_STATION_LAT", "0.0")) * u.deg
STATION_LON = float(os.getenv("ROMAN_STATION_LON", "0.0")) * u.deg
STATION_ELEV = float(os.getenv("ROMAN_STATION_ELEV", "0.0")) * u.m
OBSERVER_LOC = EarthLocation(lat=STATION_LAT, lon=STATION_LON, height=STATION_ELEV)

SMTP_SERVER = os.getenv("ROMAN_SMTP_SERVER", "smtp.gmail.com")
SENDER_EMAIL = os.getenv("ROMAN_SENDER_EMAIL", "").strip()
RECIPIENT_EMAIL = os.getenv("ROMAN_RECIPIENT_EMAIL", SENDER_EMAIL).strip()

DEFAULT_PX1_X = int(os.getenv("ROMAN_PX1_X", "2044"))
DEFAULT_PX1_Y = int(os.getenv("ROMAN_PX1_Y", "2044"))
DEFAULT_PX2_X = int(os.getenv("ROMAN_PX2_X", "1024"))
DEFAULT_PX2_Y = int(os.getenv("ROMAN_PX2_Y", "1024"))

FILTER_WAVELENGTHS = {
    "F062": 0.62, "F087": 0.87, "F106": 1.06, "F129": 1.29,
    "F146": 1.46, "F158": 1.58, "F184": 1.84, "F213": 2.13,
    "F090W": 0.90, "F115W": 1.15, "F150W": 1.50, "F187N": 1.87,
    "F200W": 2.00, "F277W": 2.77, "F335M": 3.35, "F356W": 3.56,
    "F444W": 4.44, "F470N": 4.70, "F770W": 7.70, "F1130W": 11.30
}

FILTER_PLAIN_GUIDE = {
    "F062": "Red Visible Light (sees hot stars and glowing gas)",
    "F087": "Near-Infrared (pierces thin cosmic haze)",
    "F106": "Short Infrared (spots young, hot blue stars)",
    "F129": "Mid-Short Infrared (maps distant galaxy shapes)",
    "F146": "Wide Infrared (captures maximum starlight for planet hunting)",
    "F158": "Medium Infrared (cuts through thick dust clouds)",
    "F184": "Long Infrared (reveals cool stars and warm cosmic dust)",
    "F213": "Deepest Roman Infrared (peers straight into dusty star nurseries)",
    "F444W": "Deep Heat-Vision Infrared (sees glowing dust hidden from Hubble)"
}


def load_app_password():
    pwd = os.getenv("ROMAN_APP_PASS", "").strip()
    if pwd:
        return pwd
    secret_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "roman_secret.txt")
    if os.path.exists(secret_file):
        with open(secret_file, "r", encoding="utf-8") as f:
            return f.read().strip()
    return ""


ROMAN_TARGETS = [
    {
        "name": "Roman High-Latitude Wide Area Survey",
        "plain_goal": "Mapping millions of distant galaxies to uncover how Dark Energy is expanding the universe.",
        "coord": "03h32m38s -27d48m00s",
        "filter": "F184 / F158 / F106",
        "collection": "ROMAN",
        "all_sky": False
    },
    {
        "name": "Roman Galactic Bulge Time-Domain Core",
        "plain_goal": "Watching hundreds of millions of stars at the Milky Way's center to catch hidden exoplanets and rogue worlds.",
        "coord": "17h53m24s -29d49m32s",
        "filter": "F213 / F146 / F087",
        "collection": "ROMAN",
        "all_sky": False
    },
    {
        "name": "LIVE CALIBRATION: Carina Nebula (JWST MAST)",
        "plain_goal": "Peering through thick clouds of star-forming gas and dust using real NASA infrared space telescope data.",
        "coord": "10h43m57s -59d33m51s",
        "filter": "F444W / F200W / F090W (REAL FITS)",
        "collection": "JWST",
        "all_sky": False
    },
    {
        "name": "Latest Public L2 Release (Roman All-Sky)",
        "plain_goal": "Live first-look check of the newest public deep-space image downlinked from NASA's Roman Space Telescope.",
        "coord": "03h32m38s -27d48m00s",
        "filter": "AUTO RGB (MAST Live)",
        "collection": "ROMAN",
        "all_sky": True
    }
]

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
CACHE_DIR = os.path.join(BASE_DIR, "temp_cache_live")
REPORTS_DIR = os.path.join(BASE_DIR, "reports_live")
os.makedirs(CACHE_DIR, exist_ok=True)
os.makedirs(REPORTS_DIR, exist_ok=True)


def translate_science_to_plain_english(snr, peak_flux, median_sky, px1_val, px2_val, filter_band, eph, is_live):
    """Converts raw detector numbers and sky coordinates into simple, human-readable insights."""
    if snr >= 100:
        clarity_msg = f"Crystal Clear ({snr:.1f}x brighter than background noise - ultra-sharp scientific detail)"
    elif snr >= 30:
        clarity_msg = f"Very Clean ({snr:.1f}x above background static - stars and dust clouds stand out clearly)"
    elif snr >= 10:
        clarity_msg = f"Good Signal ({snr:.1f}x above background static - faint cosmic structures visible)"
    else:
        clarity_msg = f"Faint / Noisy ({snr:.1f}x above static - capturing very dim deep-space light)"

    def describe_pixel(val):
        ratio = val / max(median_sky, 0.1)
        if ratio >= 5.0:
            return f"Direct hit on a bright star or dense glowing gas core ({ratio:.1f}x brighter than empty space)"
        elif ratio >= 1.5:
            return f"Inside a glowing nebula cloud or faint galaxy halo ({ratio:.1f}x brighter than empty space)"
        else:
            return "Clear patch of dark interstellar space (ideal baseline for measuring background light)"

    px1_plain = describe_pixel(px1_val)
    px2_plain = describe_pixel(px2_val)

    filter_desc = "Infrared heat-vision sensors (revealing stars and dust invisible to human eyes)"
    for code, desc in FILTER_PLAIN_GUIDE.items():
        if code in str(filter_band).upper():
            filter_desc = desc
            break

    if eph:
        if eph["alt"] > 0:
            sky_msg = f"Currently above your local horizon ({eph['alt']:.1f} deg elevation, azimuth {eph['az']:.0f} deg)."
        else:
            sky_msg = f"Currently on the opposite side of Earth ({abs(eph['alt']):.1f} deg below local horizon), while Roman observes freely from deep space."
    else:
        sky_msg = "Deep-space target coordinates locked."

    source_plain = (
        "Real space telescope data downloaded live from NASA's MAST archive."
        if is_live
        else "Pre-launch full-color test simulation (ready for Roman's Dec 7 live data drop)."
    )

    return {
        "source_plain": source_plain,
        "clarity_msg": clarity_msg,
        "px1_plain": px1_plain,
        "px2_plain": px2_plain,
        "filter_desc": filter_desc,
        "sky_msg": sky_msg
    }


def compute_ephemeris(coord_str=None, ra_deg=None, dec_deg=None):
    now_utc = Time.now()
    try:
        if ra_deg is not None and dec_deg is not None:
            coord = SkyCoord(ra=float(ra_deg) * u.deg, dec=float(dec_deg) * u.deg)
        elif coord_str:
            coord = SkyCoord(coord_str)
        else:
            return None

        altaz = coord.transform_to(AltAz(obstime=now_utc, location=OBSERVER_LOC))
        alt = altaz.alt.degree
        az = altaz.az.degree
        airmass = altaz.secz.value if alt > 0 else None
        status = "ABOVE HORIZON [VISIBLE FROM STATION]" if alt > 0 else "BELOW HORIZON [OTHER SIDE OF EARTH]"
        return {
            "ra": coord.ra.to_string(unit=u.hour, precision=1),
            "dec": coord.dec.to_string(unit=u.deg, precision=1),
            "alt": alt,
            "az": az,
            "airmass": airmass,
            "status": status,
            "time": now_utc.iso
        }
    except Exception:
        return None


def generate_synthetic_wfi_rgb(target_name):
    np.random.seed(abs(hash(target_name)) % (2**32))
    size = 1024
    y_grid, x_grid = np.ogrid[:size, :size]

    ch_r = np.random.poisson(lam=115, size=(size, size)).astype(np.float32)
    ch_g = np.random.poisson(lam=120, size=(size, size)).astype(np.float32)
    ch_b = np.random.poisson(lam=110, size=(size, size)).astype(np.float32)

    cx1, cy1 = size * 0.42, size * 0.54
    dist1 = np.sqrt(((x_grid - cx1) * 0.85) ** 2 + ((y_grid - cy1) * 1.15) ** 2)
    nebula_core = 260.0 * np.exp(-(dist1 ** 2) / (2 * (175.0 ** 2)))

    cx2, cy2 = size * 0.62, size * 0.38
    dist2 = np.sqrt(((x_grid - cx2) * 1.2) ** 2 + ((y_grid - cy2) * 0.8) ** 2)
    nebula_ridge = 190.0 * np.exp(-(dist2 ** 2) / (2 * (135.0 ** 2)))

    wave_mod = 1.0 + 0.25 * np.sin(x_grid / 55.0 + y_grid / 70.0)
    ch_r += (nebula_core * 1.35 + nebula_ridge * 0.85) * wave_mod
    ch_g += (nebula_core * 0.90 + nebula_ridge * 1.10) * wave_mod
    ch_b += (nebula_core * 0.45 + nebula_ridge * 1.30) * wave_mod

    num_sources = 280
    ys = np.random.randint(25, size - 25, num_sources)
    xs = np.random.randint(25, size - 25, num_sources)
    fluxes = np.random.exponential(scale=650, size=num_sources) + 90
    spectral_types = np.random.choice(["hot_blue", "solar_gold", "red_giant"], size=num_sources, p=[0.30, 0.40, 0.30])

    for y, x, f, stype in zip(ys, xs, fluxes, spectral_types):
        rad = np.random.choice([1, 2, 3], p=[0.55, 0.30, 0.15])
        y1, y2 = max(0, y - rad), min(size, y + rad + 1)
        x1, x2 = max(0, x - rad), min(size, x + rad + 1)
        amp = f / (rad + 1)

        if stype == "hot_blue":
            w_r, w_g, w_b = 0.55, 0.95, 1.45
        elif stype == "solar_gold":
            w_r, w_g, w_b = 1.15, 1.05, 0.70
        else:
            w_r, w_g, w_b = 1.55, 0.85, 0.45

        ch_r[y1:y2, x1:x2] += amp * w_r
        ch_g[y1:y2, x1:x2] += amp * w_g
        ch_b[y1:y2, x1:x2] += amp * w_b

    return {"R": ch_r, "G": ch_g, "B": ch_b}


def filter_wavelength_rank(filt_name):
    if not filt_name:
        return 1.50
    upper = str(filt_name).upper()
    for code, wl in FILTER_WAVELENGTHS.items():
        if code in upper:
            return wl
    return 1.50

def fetch_mast_exposures(target_idx, max_filters=3):
    target_info = ROMAN_TARGETS[target_idx].copy()
    collection = target_info.get("collection", "ROMAN")
    is_all_sky = target_info.get("all_sky", False)
    print(f"\n[INTERCEPT] Querying {collection} pipeline stream on MAST for {target_info['name']}...", flush=True)

    downloaded_frames = []
    seen_filters = set()

    # Fast-path: reuse any already cached science files for JWST calibration
    if collection == "JWST":
        cached_fits = sorted([
            os.path.join(CACHE_DIR, f) for f in os.listdir(CACHE_DIR)
            if f.lower().endswith('_i2d.fits') and os.path.getsize(os.path.join(CACHE_DIR, f)) > 100000
        ])
        if cached_fits:
            print(f"[CACHE HIT] Loaded {len(cached_fits[:max_filters])} cached JWST frames instantly.", flush=True)
            return [(p, "") for p in cached_fits[:max_filters]], target_info

    try:
        from astroquery.mast import Observations
        Observations.TIMEOUT = 15

        if is_all_sky:
            obs = Observations.query_criteria(
                obs_collection=collection,
                dataproduct_type="image",
                t_min=[Time('2026-11-15').mjd, Time.now().mjd + 1]
            )
            if len(obs) > 0 and 't_min' in obs.colnames:
                obs.sort('t_min')
                obs = obs[::-1]
        elif collection == "JWST":
            coord = SkyCoord(target_info["coord"])
            obs = Observations.query_criteria(
                obs_collection="JWST",
                coordinates=coord,
                radius=0.05 * u.deg,
                dataproduct_type="image",
                calib_level=3
            )
        else:
            coord = SkyCoord(target_info["coord"])
            obs = Observations.query_criteria(
                obs_collection=collection,
                coordinates=coord,
                radius=0.05 * u.deg,
                dataproduct_type="image"
            )

        if len(obs) > 0:
            for row in obs[:10]:
                if len(downloaded_frames) >= max_filters:
                    break

                obsid = str(row['obsid'])
                row_filt = str(row['filters']).upper() if 'filters' in obs.colnames and row['filters'] else ""
                if row_filt and row_filt in seen_filters and len(obs) > max_filters:
                    continue

                target_label = str(row['target_name']) if 'target_name' in obs.colnames and row['target_name'] else target_info['name']
                rest_url = (
                    f"https://mast.stsci.edu/api/v0/invoke?request="
                    f"%7B%22service%22%3A%22Mast.Caom.Products%22%2C%22params%22%3A"
                    f"%7B%22obsid%22%3A%22{obsid}%22%7D%2C%22format%22%3A%22json%22%7D"
                )
                req = urllib.request.Request(rest_url, headers={'User-Agent': 'Mozilla/5.0'})
                with urllib.request.urlopen(req, timeout=12) as resp:
                    p_data = json.loads(resp.read().decode('utf-8'))

                for prod in p_data.get('data', []):
                    fname = str(prod.get('productFilename', '')).lower()
                    ptype = str(prod.get('productType', '')).upper()
                    is_valid_ext = fname.endswith(('.asdf', '.fits', '.fits.gz'))
                    if collection == "JWST":
                        is_valid_ext = fname.endswith('_i2d.fits')
                    if ptype == "SCIENCE" and is_valid_ext:
                        safe_fname = os.path.basename(prod['productFilename'])
                        out_path = os.path.join(CACHE_DIR, safe_fname)
                        if not os.path.exists(out_path) or os.path.getsize(out_path) < 100000:
                            stream_url = f"https://mast.stsci.edu/api/v0.1/Download/file?uri={prod['dataURI']}"
                            s_req = urllib.request.Request(stream_url, headers={'User-Agent': 'Mozilla/5.0'})
                            with urllib.request.urlopen(s_req, timeout=20) as s_resp, open(out_path, 'wb') as out_f:
                                out_f.write(s_resp.read())
                        if os.path.exists(out_path) and os.path.getsize(out_path) > 100000:
                            if is_all_sky and target_label:
                                target_info["name"] = f"L2 Live: {target_label}"
                            downloaded_frames.append((out_path, row_filt))
                            if row_filt:
                                seen_filters.add(row_filt)
                            break
    except Exception as e:
        print(f"[DOWNLINK STATUS] MAST live query timed out/deferred ({e}) -> switching to fast local render.", flush=True)

    return downloaded_frames, target_info


def extract_asdf_metadata(meta_tree):
    meta_dict = {}
    if meta_tree is None:
        return meta_dict
    try:
        inst = getattr(meta_tree, 'instrument', None)
        if inst is not None:
            opt_elem = getattr(inst, 'optical_element', None)
            if opt_elem:
                meta_dict['filter'] = str(opt_elem)

        exp = getattr(meta_tree, 'exposure', None)
        if exp is not None:
            eff_exp = getattr(exp, 'effective_exposure_time', None) or getattr(exp, 'exposure_time', None)
            if eff_exp is not None:
                meta_dict['exptime'] = f"{float(eff_exp):.2f}"

        obs = getattr(meta_tree, 'observation', None)
        if obs is not None:
            start_t = getattr(obs, 'start_time', None)
            if start_t is not None:
                meta_dict['date_obs'] = str(start_t).replace('T', ' ')[:19]

        tgt = getattr(meta_tree, 'target', None)
        if tgt is not None:
            ra = getattr(tgt, 'ra', None)
            dec = getattr(tgt, 'dec', None)
            if ra is not None and dec is not None:
                meta_dict['ra_deg'] = float(ra)
                meta_dict['dec_deg'] = float(dec)
            tname = getattr(tgt, 'proposer_name', None) or getattr(tgt, 'catalog_name', None)
            if tname:
                meta_dict['target_name'] = str(tname)
    except Exception:
        pass
    return meta_dict


def extract_fits_metadata(header):
    meta_dict = {}
    if header is None:
        return meta_dict
    try:
        filt = header.get('FILTER') or header.get('FILTER1') or header.get('OPT_ELEM')
        if filt:
            meta_dict['filter'] = str(filt)
        exptime = header.get('EXPTIME') or header.get('EFFEXPTM')
        if exptime is not None:
            meta_dict['exptime'] = f"{float(exptime):.2f}"
        date_obs = header.get('DATE-OBS')
        if date_obs:
            meta_dict['date_obs'] = str(date_obs).replace('T', ' ')[:19]
        ra_val = header.get('RA_TARG') or header.get('CRVAL1')
        dec_val = header.get('DEC_TARG') or header.get('CRVAL2')
        if ra_val is not None and dec_val is not None:
            meta_dict['ra_deg'] = float(ra_val)
            meta_dict['dec_deg'] = float(dec_val)
        targ = header.get('TARGNAME') or header.get('OBJECT')
        if targ:
            meta_dict['target_name'] = str(targ)
    except Exception:
        pass
    return meta_dict


def extract_binned_array(file_path, max_dim=1400):
    if not file_path or not os.path.exists(file_path):
        return None, {}

    if file_path.lower().endswith('.asdf'):
        try:
            import roman_datamodels as rdm
            with rdm.open(file_path) as dm:
                if hasattr(dm, 'data') and dm.data is not None:
                    raw_arr = np.asarray(dm.data)
                    if raw_arr.ndim == 2:
                        step = max(1, max(raw_arr.shape) // max_dim)
                        binned = np.array(raw_arr[::step, ::step], dtype=np.float32)
                    elif raw_arr.ndim == 3 and raw_arr.shape[0] > 0:
                        step = max(1, max(raw_arr.shape[1:]) // max_dim)
                        binned = np.array(raw_arr[-1, ::step, ::step], dtype=np.float32)
                    else:
                        return None, {}
                    meta_dict = extract_asdf_metadata(getattr(dm, 'meta', None))
                    return binned, meta_dict
        except Exception:
            pass
        return None, {}

    try:
        with fits.open(file_path, memmap=True) as hdul:
            primary_meta = extract_fits_metadata(hdul[0].header) if len(hdul) > 0 else {}
            for hdu in hdul:
                if hdu.data is not None and isinstance(hdu.data, np.ndarray) and hdu.data.size > 10000:
                    meta_dict = primary_meta.copy()
                    meta_dict.update(extract_fits_metadata(hdu.header))
                    if hdu.data.ndim == 2:
                        step = max(1, max(hdu.data.shape) // max_dim)
                        return np.array(hdu.data[::step, ::step], dtype=np.float32), meta_dict
                    elif hdu.data.ndim == 3 and hdu.data.shape[0] > 0:
                        step = max(1, max(hdu.data.shape[1:]) // max_dim)
                        return np.array(hdu.data[-1, ::step, ::step], dtype=np.float32), meta_dict
    except Exception:
        pass
    return None, {}


def resize_to_match(arr, target_shape):
    if arr.shape == target_shape:
        return arr
    pil_f = Image.fromarray(arr.astype(np.float32), mode='F')
    pil_r = pil_f.resize((target_shape[1], target_shape[0]), Image.Resampling.BILINEAR)
    return np.asarray(pil_r, dtype=np.float32)


def build_lupton_composite(ch_r, ch_g, ch_b):
    def prep_channel(arr):
        clean = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
        positive = clean[clean > 0]
        sky = float(np.median(positive)) if len(positive) > 0 else 0.0
        sub = np.clip(clean - sky, 0.0, None)
        p_high = float(np.percentile(sub, 99.6)) if np.max(sub) > 0 else 1.0
        return sub / max(p_high, 1e-5)

    r_norm = prep_channel(ch_r)
    g_norm = prep_channel(ch_g)
    b_norm = prep_channel(ch_b)

    return make_lupton_rgb(r_norm, g_norm, b_norm, minimum=0.0, stretch=0.45, Q=8)


def dispatch_email(target, plain_goal, filter_band, exptime, date_obs, eph, peak_flux, median_sky, snr,
                   px1_coord, px1_val, px2_coord, px2_val, plain_summary, img_path, is_live, color_mode, collection):
    app_password = load_app_password()
    if not SENDER_EMAIL or not RECIPIENT_EMAIL or not app_password:
        return False, "Report saved locally (set ROMAN_SENDER_EMAIL & ROMAN_APP_PASS to enable email)"

    body = f"""======================================================================
     NANCY GRACE ROMAN SPACE TELESCOPE - SIMPLE SCIENCE BRIEFING
======================================================================
TARGET           : {target}
DATA STATUS      : {plain_summary['source_plain']}
MISSION GOAL     : {plain_goal}

1. WHAT THIS IMAGE SHOWS (IN SIMPLE WORDS)
----------------------------------------------------------------------
- Camera Vision    : {plain_summary['filter_desc']}
- Picture Quality  : {plain_summary['clarity_msg']}
- Exposure Length  : Captured over {exptime} seconds on {date_obs} UTC
- Sky Location     : {plain_summary['sky_msg']}

2. TRACKED DETECTOR PIXELS (WHAT THEY ARE LOOKING AT)
----------------------------------------------------------------------
- Pixel #1 (Cyan Reticle at X={px1_coord[0]}, Y={px1_coord[1]}):
  Brightness = {px1_val:.2f} light units -> {plain_summary['px1_plain']}

- Pixel #2 (Gold Reticle at X={px2_coord[0]}, Y={px2_coord[1]}):
  Brightness = {px2_val:.2f} light units -> {plain_summary['px2_plain']}

3. QUICK BRIGHTNESS COMPARISON
----------------------------------------------------------------------
- Brightest Star / Core : {peak_flux:.1f} light units
- Empty Dark Space      : {median_sky:.2f} light units
- Sky Coordinates       : RA {eph['ra'] if eph else 'N/A'} | DEC {eph['dec'] if eph else 'N/A'}
- Render Pipeline       : {color_mode} ({filter_band})
======================================================================
"""
    msg = MIMEMultipart()
    msg['From'] = SENDER_EMAIL
    msg['To'] = RECIPIENT_EMAIL
    msg['Date'] = formatdate(localtime=True)
    msg['Message-ID'] = make_msgid()
    msg['Subject'] = f"[ROMAN BRIEFING] {target} - Simple Science Report"
    msg.attach(MIMEText(body, 'plain'))

    with open(img_path, 'rb') as f:
        mime_img = MIMEImage(f.read())
        mime_img.add_header('Content-Disposition', 'attachment', filename=os.path.basename(img_path))
        msg.attach(mime_img)

    try:
        with smtplib.SMTP(SMTP_SERVER, 587, timeout=15) as server:
            server.ehlo()
            server.starttls()
            server.ehlo()
            server.login(SENDER_EMAIL, app_password)
            server.send_message(msg)
        return True, "Delivered via Port 587"
    except Exception:
        try:
            with smtplib.SMTP_SSL(SMTP_SERVER, 465, timeout=15) as server:
                server.login(SENDER_EMAIL, app_password)
                server.send_message(msg)
            return True, "Delivered via Port 465"
        except Exception as e465:
            return False, str(e465)


# =====================================================================
# WORKSTATION INTERFACE
# =====================================================================
class NancyGroundCommand(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("NANCY GRACE ROMAN SPACE TELESCOPE // GROUND TELEMETRY WORKSTATION [v3 BETA 2]")
        self.geometry("1260x850")
        self.configure(bg="#04070d")

        self.create_widgets()
        self.draw_radar(0, 0, "STANDBY")
        self.update_clock()

    def create_widgets(self):
        topbar = tk.Frame(self, bg="#080e18", height=54, bd=1, relief="solid")
        topbar.pack(fill="x", padx=10, pady=(8, 4))

        tk.Label(topbar, text="NANCY GRACE ROMAN [v3 BETA 2]", font=("Consolas", 15, "bold"), fg="#00e5ff", bg="#080e18").pack(side="left", padx=15)
        tk.Label(topbar, text="| 300-MEGAPIXEL WIDE FIELD GROUND STATION // SIMPLE SCIENCE TRANSLATOR", font=("Consolas", 9), fg="#64748b", bg="#080e18").pack(side="left")

        self.time_lbl = tk.Label(topbar, text="", font=("Consolas", 10, "bold"), fg="#10b981", bg="#080e18")
        self.time_lbl.pack(side="right", padx=15)

        deck = tk.Frame(self, bg="#04070d")
        deck.pack(fill="both", expand=True, padx=10, pady=4)

        left_pnl = tk.Frame(deck, bg="#080e18", width=365, bd=1, relief="solid")
        left_pnl.pack(side="left", fill="y", padx=(0, 6), pady=4)

        tk.Label(left_pnl, text="CHOOSE SPACE TARGET", font=("Consolas", 10, "bold"), fg="#f59e0b", bg="#080e18").pack(pady=(10, 6))

        self.btn_target0 = tk.Button(left_pnl, text="[1] Deep Galaxy & Dark Energy Map", font=("Consolas", 9, "bold"),
                                     fg="#00e5ff", bg="#101927", activebackground="#00e5ff", relief="flat", padx=10, pady=6,
                                     command=lambda: self.launch_acquisition(0))
        self.btn_target0.pack(fill="x", padx=12, pady=3)

        self.btn_target1 = tk.Button(left_pnl, text="[2] Milky Way Core (Planet Hunter)", font=("Consolas", 9, "bold"),
                                     fg="#38bdf8", bg="#101927", activebackground="#38bdf8", relief="flat", padx=10, pady=6,
                                     command=lambda: self.launch_acquisition(1))
        self.btn_target1.pack(fill="x", padx=12, pady=3)

        self.btn_target2 = tk.Button(left_pnl, text="[3] LIVE TEST: Carina Star Nursery", font=("Consolas", 9, "bold"),
                                     fg="#f43f5e", bg="#1f1218", activebackground="#f43f5e", relief="flat", padx=10, pady=6,
                                     command=lambda: self.launch_acquisition(2))
        self.btn_target2.pack(fill="x", padx=12, pady=3)

        self.btn_target3 = tk.Button(left_pnl, text="[4] LATEST LIVE ROMAN RELEASE", font=("Consolas", 9, "bold"),
                                     fg="#fbbf24", bg="#1c1917", activebackground="#fbbf24", relief="flat", padx=10, pady=6,
                                     command=lambda: self.launch_acquisition(3))
        self.btn_target3.pack(fill="x", padx=12, pady=3)

        pix_box = tk.LabelFrame(left_pnl, text=" DUAL PIXEL RETICLE TARGETER ", font=("Consolas", 8, "bold"),
                                fg="#00e5ff", bg="#080e18", bd=1, relief="solid")
        pix_box.pack(fill="x", padx=12, pady=(10, 6))

        row1 = tk.Frame(pix_box, bg="#080e18")
        row1.pack(fill="x", padx=6, pady=(4, 2))
        tk.Label(row1, text="PX1 X:", font=("Consolas", 8, "bold"), fg="#00e5ff", bg="#080e18").pack(side="left")
        self.entry_x1 = tk.Entry(row1, font=("Consolas", 8), width=7, bg="#04070d", fg="#00e5ff", insertbackground="#00e5ff")
        self.entry_x1.insert(0, str(DEFAULT_PX1_X))
        self.entry_x1.pack(side="left", padx=4)

        tk.Label(row1, text="Y:", font=("Consolas", 8, "bold"), fg="#00e5ff", bg="#080e18").pack(side="left")
        self.entry_y1 = tk.Entry(row1, font=("Consolas", 8), width=7, bg="#04070d", fg="#00e5ff", insertbackground="#00e5ff")
        self.entry_y1.insert(0, str(DEFAULT_PX1_Y))
        self.entry_y1.pack(side="left", padx=4)

        row2 = tk.Frame(pix_box, bg="#080e18")
        row2.pack(fill="x", padx=6, pady=(2, 4))
        tk.Label(row2, text="PX2 X:", font=("Consolas", 8, "bold"), fg="#fbbf24", bg="#080e18").pack(side="left")
        self.entry_x2 = tk.Entry(row2, font=("Consolas", 8), width=7, bg="#04070d", fg="#fbbf24", insertbackground="#fbbf24")
        self.entry_x2.insert(0, str(DEFAULT_PX2_X))
        self.entry_x2.pack(side="left", padx=4)

        tk.Label(row2, text="Y:", font=("Consolas", 8, "bold"), fg="#fbbf24", bg="#080e18").pack(side="left")
        self.entry_y2 = tk.Entry(row2, font=("Consolas", 8), width=7, bg="#04070d", fg="#fbbf24", insertbackground="#fbbf24")
        self.entry_y2.insert(0, str(DEFAULT_PX2_Y))
        self.entry_y2.pack(side="left", padx=4)

        self.lbl_pix1_status = tk.Label(pix_box, text=f"PX1 ({DEFAULT_PX1_X}, {DEFAULT_PX1_Y}): Ready", font=("Consolas", 7),
                                        fg="#00e5ff", bg="#080e18")
        self.lbl_pix1_status.pack(pady=(1, 1))

        self.lbl_pix2_status = tk.Label(pix_box, text=f"PX2 ({DEFAULT_PX2_X}, {DEFAULT_PX2_Y}): Ready", font=("Consolas", 7),
                                        fg="#fbbf24", bg="#080e18")
        self.lbl_pix2_status.pack(pady=(1, 5))

        tk.Label(left_pnl, text="LOCAL SKY COMPASS", font=("Consolas", 9, "bold"), fg="#10b981", bg="#080e18").pack(pady=(6, 2))
        self.radar_canvas = tk.Canvas(left_pnl, width=175, height=175, bg="#04070d", highlightthickness=1, highlightbackground="#1e293b")
        self.radar_canvas.pack(pady=4)

        self.lbl_pos = tk.Label(left_pnl, text="SKY COORD: --:--:-- | --°--'--", font=("Consolas", 8), fg="#64748b", bg="#080e18")
        self.lbl_pos.pack(pady=(4, 1))
        self.lbl_altaz = tk.Label(left_pnl, text="HEIGHT: --.-° | COMPASS: ---.-°", font=("Consolas", 8), fg="#64748b", bg="#080e18")
        self.lbl_altaz.pack(pady=1)
        self.lbl_snr = tk.Label(left_pnl, text="CLARITY: --.--x | PEAK: --", font=("Consolas", 8, "bold"), fg="#00e5ff", bg="#080e18")
        self.lbl_snr.pack(pady=1)

        right_pnl = tk.Frame(deck, bg="#080e18", bd=1, relief="solid")
        right_pnl.pack(side="right", fill="both", expand=True, padx=(6, 0), pady=4)

        status_strip = tk.Frame(right_pnl, bg="#0d1524", height=38)
        status_strip.pack(fill="x", padx=8, pady=(8, 4))

        self.status_msg = tk.StringVar(value="WORKSTATION READY // SELECT A SPACE TARGET ON THE LEFT")
        self.status_lbl = tk.Label(status_strip, textvariable=self.status_msg, font=("Consolas", 9, "bold"), fg="#10b981", bg="#0d1524")
        self.status_lbl.pack(side="left", padx=10)

        self.image_display = tk.Label(right_pnl, bg="#000000",
                                      text="SPACE TELESCOPE VIEWPORT\nClick any button on the left to download and translate deep-space data",
                                      font=("Consolas", 11), fg="#475569")
        self.image_display.pack(fill="both", expand=True, padx=8, pady=4)

        brief_box = tk.LabelFrame(right_pnl, text=" PLAIN-ENGLISH SCIENCE TRANSLATION ", font=("Consolas", 9, "bold"),
                                  fg="#10b981", bg="#0d1524", bd=1, relief="solid")
        brief_box.pack(fill="x", padx=8, pady=(4, 8))

        self.brief_text = tk.StringVar(
            value="Waiting for observation... Once data loads, this box explains in simple words what the telescope sees, how clear the image is, and what your 2 pixels hit."
        )
        self.brief_lbl = tk.Label(brief_box, textvariable=self.brief_text, font=("Consolas", 8), fg="#e2e8f0",
                                  bg="#0d1524", justify="left", anchor="w", wraplength=820)
        self.brief_lbl.pack(fill="x", padx=10, pady=8)

    def update_clock(self):
        utc_time = datetime.now(timezone.utc).strftime('%Y-%m-%d  %H:%M:%S UTC')
        self.time_lbl.config(text=utc_time)
        self.after(1000, self.update_clock)

    def set_status(self, text, color="#10b981"):
        self.after(0, lambda: (self.status_msg.set(text), self.status_lbl.config(fg=color)))

    def draw_radar(self, az, alt, status):
        self.radar_canvas.delete("all")
        cx, cy, r = 87, 87, 70

        for ring in [r, int(r * 0.66), int(r * 0.33)]:
            self.radar_canvas.create_oval(cx - ring, cy - ring, cx + ring, cy + ring, outline="#1e293b", width=1)
        self.radar_canvas.create_line(cx - r, cy, cx + r, cy, fill="#1e293b")
        self.radar_canvas.create_line(cx, cy - r, cx, cy + r, fill="#1e293b")

        self.radar_canvas.create_text(cx, cy - r + 8, text="N", fill="#00e5ff", font=("Consolas", 7, "bold"))
        self.radar_canvas.create_text(cx + r - 8, cy, text="E", fill="#00e5ff", font=("Consolas", 7, "bold"))
        self.radar_canvas.create_text(cx, cy + r - 8, text="S", fill="#00e5ff", font=("Consolas", 7, "bold"))
        self.radar_canvas.create_text(cx - r + 8, cy, text="W", fill="#00e5ff", font=("Consolas", 7, "bold"))

        if status != "STANDBY":
            rad = math.radians(az - 90)
            dist = r * (1.0 - max(0, alt) / 90.0)
            px = cx + dist * math.cos(rad)
            py = cy + dist * math.sin(rad)
            color = "#10b981" if alt > 0 else "#f43f5e"

            self.radar_canvas.create_line(cx, cy, px, py, fill=color, width=2)
            self.radar_canvas.create_oval(px - 5, py - 5, px + 5, py + 5, fill=color, outline="#ffffff")

    def launch_acquisition(self, idx):
        self.btn_target0.config(state="disabled")
        self.btn_target1.config(state="disabled")
        self.btn_target2.config(state="disabled")
        self.btn_target3.config(state="disabled")
        threading.Thread(target=self.run_pipeline, args=(idx,), daemon=True).start()

    def run_pipeline(self, idx):
        target_info = ROMAN_TARGETS[idx]
        collection = target_info.get("collection", "ROMAN")
        plain_goal = target_info.get("plain_goal", "Deep-space infrared observation.")
        self.set_status(f"CONNECTING TO NASA ARCHIVE FOR {target_info['name']}...", "#38bdf8")

        eph = compute_ephemeris(coord_str=target_info["coord"])
        if eph:
            self.after(0, lambda: self.update_telemetry_hud(eph))

        downloaded_frames, meta_target = fetch_mast_exposures(idx, max_filters=3)

        parsed_bands = []
        for fpath, row_filt in downloaded_frames:
            arr, hdr_meta = extract_binned_array(fpath, max_dim=1200)
            if arr is not None:
                filt_label = hdr_meta.get("filter") or row_filt or "WFI-IR"
                parsed_bands.append({
                    "array": arr,
                    "filter": filt_label,
                    "wl": filter_wavelength_rank(filt_label),
                    "meta": hdr_meta
                })

        is_live = len(parsed_bands) > 0

        if not is_live:
            self.set_status("BUILDING 3-COLOR TEST VIEW (WAITING FOR ROMAN DEC 7 RELEASE)...", "#38bdf8")
            rgb_dict = generate_synthetic_wfi_rgb(meta_target["name"])
            ch_r, ch_g, ch_b = rgb_dict["R"], rgb_dict["G"], rgb_dict["B"]
            lum_data = (ch_r + ch_g + ch_b) / 3.0
            display_target = meta_target["name"]
            filter_band = meta_target.get("filter", "F184 / F158 / F106")
            exptime = "142.5"
            date_obs = datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
            color_mode = "3-Color Infrared Blend (Test Mode)"
            rgb_img = build_lupton_composite(ch_r, ch_g, ch_b)
            is_rgb = True
        else:
            parsed_bands.sort(key=lambda x: x["wl"])
            primary = parsed_bands[-1]
            hdr_meta = primary["meta"]

            display_target = hdr_meta.get("target_name") or meta_target["name"]
            exptime = hdr_meta.get("exptime", "142.5")
            date_obs = hdr_meta.get("date_obs") or datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')

            if "ra_deg" in hdr_meta and "dec_deg" in hdr_meta:
                live_eph = compute_ephemeris(ra_deg=hdr_meta["ra_deg"], dec_deg=hdr_meta["dec_deg"])
                if live_eph:
                    eph = live_eph
                    self.after(0, lambda: self.update_telemetry_hud(eph))

            base_shape = primary["array"].shape
            for b in parsed_bands:
                if b["array"].shape != base_shape:
                    b["array"] = resize_to_match(b["array"], base_shape)

            if len(parsed_bands) >= 3:
                ch_b = parsed_bands[0]["array"]
                ch_g = parsed_bands[1]["array"]
                ch_r = parsed_bands[-1]["array"]
                filter_band = f"{parsed_bands[-1]['filter']} + {parsed_bands[1]['filter']} + {parsed_bands[0]['filter']}"
                color_mode = f"True 3-Filter Full Color ({collection} Live)"
                rgb_img = build_lupton_composite(ch_r, ch_g, ch_b)
                lum_data = (ch_r + ch_g + ch_b) / 3.0
                is_rgb = True
            elif len(parsed_bands) == 2:
                ch_b = parsed_bands[0]["array"]
                ch_r = parsed_bands[1]["array"]
                ch_g = 0.5 * (ch_r + ch_b)
                filter_band = f"{parsed_bands[1]['filter']} + {parsed_bands[0]['filter']}"
                color_mode = f"2-Filter Color Blend ({collection} Live)"
                rgb_img = build_lupton_composite(ch_r, ch_g, ch_b)
                lum_data = 0.5 * (ch_r + ch_b)
                is_rgb = True
            else:
                lum_data = primary["array"]
                filter_band = primary["filter"]
                color_mode = f"Single-Filter Heat Map ({collection} Live)"
                rgb_img = None
                is_rgb = False

        valid = lum_data[np.isfinite(lum_data) & (lum_data > 0)]
        peak_flux = float(np.max(valid)) if len(valid) > 0 else 0.0
        median_sky = float(np.median(valid)) if len(valid) > 0 else 0.0
        snr = (peak_flux - median_sky) / np.sqrt(peak_flux) if peak_flux > 0 else 0.0

        scale_x = 4088.0 / lum_data.shape[1]
        scale_y = 4088.0 / lum_data.shape[0]

        def resolve_pixel(entry_x, entry_y, def_x, def_y):
            try:
                rx = int(entry_x.get())
                ry = int(entry_y.get())
            except Exception:
                rx, ry = def_x, def_y
            px = max(0, min(int(rx / scale_x), lum_data.shape[1] - 1))
            py = max(0, min(int(ry / scale_y), lum_data.shape[0] - 1))
            val = float(lum_data[py, px])
            return rx, ry, px, py, val

        raw_x1, raw_y1, target_px1, target_py1, pixel_val1 = resolve_pixel(self.entry_x1, self.entry_y1, DEFAULT_PX1_X, DEFAULT_PX1_Y)
        raw_x2, raw_y2, target_px2, target_py2, pixel_val2 = resolve_pixel(self.entry_x2, self.entry_y2, DEFAULT_PX2_X, DEFAULT_PX2_Y)

        plain_summary = translate_science_to_plain_english(
            snr, peak_flux, median_sky, pixel_val1, pixel_val2, filter_band, eph, is_live
        )

        self.after(0, lambda: self.lbl_pix1_status.config(
            text=f"PX1 ({raw_x1}, {raw_y1}): {pixel_val1:.1f} light units"
        ))
        self.after(0, lambda: self.lbl_pix2_status.config(
            text=f"PX2 ({raw_x2}, {raw_y2}): {pixel_val2:.1f} light units"
        ))
        self.after(0, lambda: self.lbl_snr.config(text=f"CLARITY: {snr:.1f}x | PEAK: {peak_flux:.0f}"))

        briefing_display = (
            f"WHAT THIS IS: {plain_goal}\n"
            f"CAMERA VISION: {plain_summary['filter_desc']} | CLARITY: {plain_summary['clarity_msg']}\n"
            f"PIXEL #1 ({raw_x1}, {raw_y1}): {plain_summary['px1_plain']}\n"
            f"PIXEL #2 ({raw_x2}, {raw_y2}): {plain_summary['px2_plain']}\n"
            f"SKY POSITION: {plain_summary['sky_msg']}"
        )
        self.after(0, lambda: self.brief_text.set(briefing_display))

        safe_slug = "".join(c if c.isalnum() or c in "_-" else "_" for c in display_target)
        out_img = os.path.join(REPORTS_DIR, f"ROMAN_V3_RGB_{safe_slug}.png")
        fig, (ax_img, ax_hist) = plt.subplots(1, 2, figsize=(11.5, 4.4), facecolor="#04070d")

        mode_tag = f"LIVE {collection}" if is_live else "PRE-LAUNCH TEST"
        if is_rgb and rgb_img is not None:
            ax_img.imshow(rgb_img, origin="lower")
        else:
            if len(valid) > 0:
                interval = ZScaleInterval(contrast=0.25)
                vmin, vmax = interval.get_limits(valid)
                norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(a=0.1))
            else:
                norm = None
            ax_img.imshow(lum_data, cmap="inferno", norm=norm, origin="lower")

        ax_img.set_title(f"{display_target}\n[{mode_tag} | Filter: {filter_band}]",
                         color="#00e5ff", fontsize=8.5, fontfamily="monospace")

        reticle1 = patches.Circle((target_px1, target_py1), radius=18, edgecolor="#00e5ff", facecolor="none", linewidth=1.4)
        ax_img.add_patch(reticle1)
        ax_img.plot([target_px1 - 26, target_px1 + 26], [target_py1, target_py1], color="#00e5ff", linewidth=0.9)
        ax_img.plot([target_px1, target_px1], [target_py1 - 26, target_py1 + 26], color="#00e5ff", linewidth=0.9)
        ax_img.text(target_px1 + 12, target_py1 + 12, f"PX1 ({raw_x1}, {raw_y1})",
                    color="#00e5ff", fontsize=7, fontfamily="monospace", weight="bold")

        reticle2 = patches.Circle((target_px2, target_py2), radius=18, edgecolor="#fbbf24", facecolor="none", linewidth=1.4)
        ax_img.add_patch(reticle2)
        ax_img.plot([target_px2 - 26, target_px2 + 26], [target_py2, target_py2], color="#fbbf24", linewidth=0.9)
        ax_img.plot([target_px2, target_px2], [target_py2 - 26, target_py2 + 26], color="#fbbf24", linewidth=0.9)
        ax_img.text(target_px2 + 12, target_py2 + 12, f"PX2 ({raw_x2}, {raw_y2})",
                    color="#fbbf24", fontsize=7, fontfamily="monospace", weight="bold")

        ax_img.axis("off")

        ax_hist.set_facecolor("#080e18")
        if not is_live or (is_live and is_rgb):
            for ch_arr, col, lbl in [(ch_r, "#f43f5e", "Warm Dust (Red)"), (ch_g, "#10b981", "Mid Heat (Green)"), (ch_b, "#38bdf8", "Hot Stars (Blue)")]:
                v_ch = ch_arr[np.isfinite(ch_arr) & (ch_arr > 0)]
                if len(v_ch) > 0:
                    ax_hist.hist(np.log1p(v_ch), bins=45, color=col, alpha=0.45, label=lbl, histtype="stepfilled")
            ax_hist.legend(facecolor="#04070d", edgecolor="#1e293b", labelcolor="#cbd5e1", fontsize=7)
            ax_hist.set_title("LIGHT BRIGHTNESS SPREAD (DARK -> BRIGHT)", color="#10b981", fontsize=9, fontfamily="monospace")
        else:
            if len(valid) > 0:
                ax_hist.hist(np.log1p(valid), bins=50, color="#f59e0b", edgecolor="#04070d", alpha=0.85)
            ax_hist.set_title("PIXEL BRIGHTNESS SPREAD (DARK -> BRIGHT)", color="#10b981", fontsize=9, fontfamily="monospace")
        ax_hist.set_xlabel("Brightness Level (Left = Dark Space, Right = Bright Stars)", color="#94a3b8", fontsize=7.5)
        ax_hist.set_ylabel("Number of Pixels", color="#94a3b8", fontsize=7.5)
        ax_hist.tick_params(colors="#64748b")

        plt.tight_layout()
        plt.savefig(out_img, dpi=140, facecolor=fig.get_facecolor(), edgecolor='none')
        plt.close()

        pil_img = Image.open(out_img)
        pil_img = pil_img.resize((720, 285), Image.Resampling.LANCZOS)
        tk_img = ImageTk.PhotoImage(pil_img)

        self.after(0, lambda: self.render_viewport(tk_img))

        if SENDER_EMAIL and RECIPIENT_EMAIL and load_app_password():
            self.set_status(f"SENDING SIMPLE SCIENCE REPORT TO {RECIPIENT_EMAIL}...", "#f59e0b")
        sent_ok, info = dispatch_email(
            display_target, plain_goal, filter_band, exptime, date_obs, eph,
            peak_flux, median_sky, snr,
            (raw_x1, raw_y1), pixel_val1, (raw_x2, raw_y2), pixel_val2,
            plain_summary, out_img, is_live, color_mode, collection
        )

        if sent_ok:
            self.set_status(f"DONE // SIMPLE SCIENCE REPORT EMAILED ({info})", "#10b981")
        else:
            self.set_status(f"DONE // {info}", "#10b981")

        self.after(0, self.unlock_buttons)

        del lum_data, valid
        gc.collect()

    def update_telemetry_hud(self, eph):
        self.draw_radar(eph['az'], eph['alt'], eph['status'])
        self.lbl_pos.config(text=f"SKY COORD: {eph['ra']} | {eph['dec']}")
        self.lbl_altaz.config(text=f"HEIGHT: {eph['alt']:+.1f}° | COMPASS: {eph['az']:.1f}°")

    def render_viewport(self, tk_img):
        self.image_display.configure(image=tk_img, text="")
        self.image_display.image = tk_img

    def unlock_buttons(self):
        self.btn_target0.config(state="normal")
        self.btn_target1.config(state="normal")
        self.btn_target2.config(state="normal")
        self.btn_target3.config(state="normal")


if __name__ == "__main__":
    app = NancyGroundCommand()
    app.mainloop()
