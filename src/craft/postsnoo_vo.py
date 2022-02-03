#!/usr/bin/env python
"""
Posts a snoopy log and header as a voevent to ... somewhere

Based on CREATE_VOevent by Shivani Bhandari

Copyright (C) CSIRO 2018
"""
from .crafthdr import DadaHeader
import os
import sys
import logging
import astropy.coordinates as coord
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.time import Time
import voeventparse as vp
import datetime
import pytz
import numpy as np
import argparse
from xml.dom import minidom


__author__ = "Keith Bannister <keith.bannister@csiro.au>"

obs_header = ['SOURCE','FIELD_NAME','SBID','SB_ALIAS',
              'SB_TEMPLATE','SB_TEMPLATE_VERSION',
              'SB_START_TIME','SB_OWNER','FOOTPRINT_NAME',
              'FOOTPRINT_PITCH','FOOTPRINT_ROTATION','SCANID','SCANNAME','SCAN_INTENT',
              'CAPTUREID','CAPTURE_NAME']

              

def obsparms(hdr, beamno):
    p = []
    chanbw = float(hdr.get_value('BW')) # can be negative
    nchan = int(hdr.get_value('NCHAN'))
    totalbw = abs(chanbw)*nchan
    f1 = float(hdr.get_value('FREQ'))
    centerfreq = f1 + nchan/2*chanbw
    npol = int(hdr.get_value('NPOL'))

    # parameters of the observing
    p += vp.Param(name="sampling_time", value=float(hdr.get_value('TSAMP')), unit="ms", ucd="time.resolution", ac=True)
    p += vp.Param(name="bandwidth", value=totalbw, unit="MHz", ucd="instr.bandwidth", ac=True)
    p += vp.Param(name='channel_bandwidth', value=chanbw, unit='MHz', ucd='instr.bandwidth', ac=True)
    p += vp.Param(name="nchan", value=nchan, dataType="int", ucd="meta.number;em.freq;em.bin", unit="None")
    p += vp.Param(name="centre_frequency", value=centerfreq, unit="MHz", ucd="em.freq;instr", ac=True)
    p += vp.Param(name="npol", value=npol, dataType="int", unit="None")
    p += vp.Param(name="bits_per_sample", value="8", dataType="int", unit="None")
    #p += vp.Param(name="gain", value=1.0, unit="K/Jy", ac=True) # might need to thik a bit more about his
    p += vp.Param(name="tsys", value=75.0, unit="K", ucd="phot.antennaTemp", ac=True)
    p += vp.Param(name='antennas', value=hdr.get_value('THIS_SUBARRAY'), ac=True)
    p += vp.Param(name='ics_num_ant', value=int(hdr.get_value('ICS_NUM_ANT')), ac=True)
    p += vp.Param(name="backend", value="CRAFT")

    for h in obs_header:
        p += vp.Param(name=h.lower, value=hdr.get_value(h))

    # According ot this: https://gcn.gsfc.nasa.gov/vo_examples.html should add this to wherewwhen
    ra = float(hdr.get_value('RA'))
    dec = float(hdr.get_value('DEC'))
    ant_pos = vp.Position2D(ra=ra, dec=dec, units='deg', system=vp.definitions.sky_coord_system.utc_fk5_geo)
    p += ant_pos

    beam_ra_all = list(map(float, hdr.get_value('BEAM_RA').split(',')))
    beam_dec_all = list(map(float, hdr.get_value('BEAM_DEC').split(',')))
    beam_ra = beam_ra_all[beamno]
    beam_dec = beam_dec_all[beamno]

    beam_pos = vp.Position2D(ra=beam_ra, dec=beam_dec, units='deg', system=vp.definitions.sky_coord_system.utc_fk5_geo)
    p += beam_pos
    
    group = vp.Group(params=p, name="observatory parameters")

                  
    return group


def NewVOEvent(dm, dm_err, width, snr, fluence, ra, dec,  ne2001, name, importance, utc, gl, gb): 

    z = dm/1200.0  #May change
    #errDeg = semiMaj/60.0
    error_deg = 10 # Dummy number

    # Parse UTC
    utc_YY = int(utc[:4])
    utc_MM = int(utc[5:7])
    utc_DD = int(utc[8:10])
    utc_hh = int(utc[11:13])
    utc_mm = int(utc[14:16])
    utc_ss = float(utc[17:])
    t = Time('T'.join([utc[:10], utc[11:]]), scale='utc', format='isot')
    mjd = t.mjd
    
    now = Time.now()
    mjd_now = now.mjd
   
    ivorn = ''.join([name, str(utc_hh), str(utc_mm), '/', str(mjd_now)]) 

    v = vp.Voevent(stream='aksap site address', stream_id=ivorn, role=vp.definitions.roles.test)
    #v = vp.Voevent(stream='nl.astron.apertif/alert', stream_id=ivorn, role=vp.definitions.roles.observation)
    
    # Author origin information
    #vp.set_who(v, date=datetime.datetime.utcnow(), author_ivorn="nl.astron")
    vp.set_who(v, date=datetime.datetime.utcnow(), author_ivorn="csiro?")
    # Author contact information
    vp.set_author(v, title="ASKAP FRB candidate", contactName="Shivani Bhandari", contactEmail="shivanibhandari58@gmail.com", shortName="ALERT")
    # Parameter definitions

    # These will go in the "What section"
    #ASKAP-specific observing configuration 
    #beam_sMa = vp.Param(name="beam_semi-major_axis", unit="MM", ucd="instr.beam;pos.errorEllipse;phys.angSize.smajAxis", ac=True, value=semiMaj)
    #beam_sma = vp.Param(name="beam_semi-minor_axis", unit="MM", ucd="instr.beam;pos.errorEllipse;phys.angSize.sminAxis", ac=True, value=semiMin)
    #beam_rot = vp.Param(name="beam_rotation_angle", value=0.0, unit="Degrees", ucd="instr.beam;pos.errorEllipse;instr.offset", ac=True)
    
    # Get these CRAFT observation values automatically from open_scan_rt.sum
    tsamp = vp.Param(name="sampling_time", value=0.0496, unit="ms", ucd="time.resolution", ac=True)
    bw = vp.Param(name="bandwidth", value=300.0, unit="MHz", ucd="instr.bandwidth", ac=True)
    nchan = vp.Param(name="nchan", value="1536", dataType="int", ucd="meta.number;em.freq;em.bin", unit="None")
    cf = vp.Param(name="centre_frequency", value=1400.0, unit="MHz", ucd="em.freq;instr", ac=True)
    npol = vp.Param(name="npol", value="2", dataType="int", unit="None")
    bits = vp.Param(name="bits_per_sample", value="8", dataType="int", unit="None")
    gain = vp.Param(name="gain", value=1.0, unit="K/Jy", ac=True)
    tsys = vp.Param(name="tsys", value=75.0, unit="K", ucd="phot.antennaTemp", ac=True)
    backend = vp.Param(name="backend", value="CRAFT")
#   beam = vp.Param(name="beam", value= )

    v.What.append(vp.Group(params=[tsamp, bw, nchan, cf, npol, bits, gain, tsys, backend], name="observatory parameters"))

    #Event parameters
    DM = vp.Param(name="dm", ucd="phys.dispMeasure", unit="pc/cm^3", ac=True, value=dm )
#   DM_err = vp.Param(name="dm_err", ucd="stat.error;phys.dispMeasure", unit="pc/cm^3", ac=True, value=dm_err)
    Width = vp.Param(name="width", ucd="time.duration;src.var.pulse", unit="ms", ac=True, value=width)
    SNR = vp.Param(name="snr", ucd="stat.snr", unit="None", ac=True, value=snr)
    fluence = vp.Param(name="fluence", ucd="phot.fluence", unit="Jyms", ac=True, value=fluence)
    #fluence.Description = "Calculated from radiometer equation. Not calibrated."
    Gl = vp.Param(name="gl", ucd="pos.galactic.lon", unit="Degrees", ac=True, value=gl)
    Gb = vp.Param(name="gb", ucd="pos.galactic.lat", unit="Degrees", ac=True, value=gb)

    v.What.append(vp.Group(params=[DM, Width, SNR, fluence, Gl, Gb], name="event parameters"))
#    v.What.append(vp.Group(params=[DM, DM_err, Width, SNR, fluence, Gl, Gb], name="event parameters"))

    #Advanced parameters (note, change script if using a differeing MW model)
    mw_dm = vp.Param(name="MW_dm_limit", unit="pc/cm^3", ac=True, value=ne2001)
    mw_model = vp.Param(name="galactic_electron_model", value="NE2001")
    redshift_inferred = vp.Param(name="redshift_inferred", ucd="src.redshift", unit="None", value=z)
    redshift_inferred.Description = "Redshift estimated using z = DM/1200.0 (Ioka 2003)"

    v.What.append(vp.Group(params=[mw_dm, mw_model, redshift_inferred], name="advanced parameters"))


    #WhereWhen

    vp.add_where_when(v, coords=vp.Position2D(ra=ra, dec=dec, err= error_deg, units='deg', system=vp.definitions.sky_coord_system.utc_fk5_geo),
        obs_time=datetime.datetime(utc_YY,utc_MM,utc_DD,utc_hh,utc_mm,int(utc_ss), tzinfo=pytz.UTC), observatory_location="ASKAP")

    #Why
    
    vp.add_why(v, importance=imp)
    v.Why.Name = name

    if vp.valid_as_v2_0(v):
        with open('%s.xml' % utc, 'wb') as f:
            #voxml = vp.dumps(v)
            #print (voxml)
            #xmlstr = minidom.parseString(voxml).toprettyxml(indent="   ")
            #f.write(xmlstr)
            vp.dump(v, f)
            print((vp.prettystr(v.Who)))
            print((vp.prettystr(v.What)))
            print((vp.prettystr(v.WhereWhen)))
            print((vp.prettystr(v.Why)))
    else:
        print(("Unable to write file %s.xml" % name))

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s','--snoopy-log', help='path to snoopy log', required=True, default='snoopyv2.cand')
    parser.add_argument('-d', '--dada-hdr', help='Path to Dada header')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    # S/N, sampno, secs from file start, boxcar, idt, dm, beamno, mjd, latency_ms
    #10.16 651620 1125.9994 11 998 1228.66 1 58772.0634430710 608.67
    cand = np.loadtxt(values.snoopy_log)
    sn, sampno, tsec, boxcar, idt, dm, beamno, mjd, latencyms = cand
    mjd = float(mjd)
    hdr = DadaHeader.fromfile(values.dada_hdr)

    # TODO: add fetch results to command line
    # If they are specified on the command line - partse them and add them
    # To the VO Event - otherwise don't include them in the VOEVENT

    # Move all these caluculations into the NewVoEvent function
    dm_err = 0
    width = boxcar+1
    beam_ras = list(map(float, hdr.get_value('BEAM_RA').split(',')))
    beam_decs = list(map(float, hdr.get_value('BEAM_DEC').split(',')))
    beamno = int(beamno) 
    ra = beam_ras[beamno] # degres
    dec = beam_decs[beamno] # degrees
    snr = sn
    tsamp = float(hdr.get_value('TSAMP'))
    fluence = 0 # TODO: calculate a simple fluence. Assume single antenna has an SEFD of 2000 Jy, fluence=2000/sqrt(Nant) * snr * width
    semiMag = 0 # ignore
    semiMin = 0 # ignore
    ne2001 = 0 # Find this from a python package somewher
    name= 0 # ideally we need to get this from TNS - I'm not sure what do do in the interim
    imp = '' # Is this fixed?
    utc = '2019-01-11 12:12:12' # Time of the FRB. MJD of the FRB is in the snoopy.log - it's here as the mjd variable. COnver to UTC ! BUT - this is is the time of arrival of the FRB in thebottom channel of hte band. So - we should get the frequencies from the header and work out what that channel frequency is, then correct the MJD to infinite frequency, then convert to MJD.
    # IN fact - put both times (at the bottom of the band, and at infinite frequency) and let people choose. PUt the MJDs in and UTCS in - the more the merrier.
    
    #dm = args.dm
    #dm_err = args.dm_err
    #width = args.width
    #snr = args.snr
    #fluence = args.fluence
    #ra = args.RA
    #dec = args.DEC
    #semiMaj = args.semiMaj
    #semiMin = args.semiMin
    #ne2001 = args.NE2001
    #name = args.name
    #imp = args.importance
    #utc = args.utc

    print((dm, dm_err, width, fluence, ra, dec, ne2001, name, imp, utc))

    # Parse coordinates
    c = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')
    g = c.galactic
    gl = g.l.deg
    gb = g.b.deg
   
    print((c, g, gl, gb))
    #print utc, utc_YY, utc_MM, utc_DD, utc_hh, utc_mm, utc_ss, mjd
    # TODO: put the entire cand information in from teh snoopy log.
    
    # TODO: - make this function NewVoEvent(cand, hdr)
    NewVOEvent(dm, dm_err, width, snr, fluence, ra, dec, ne2001, name, imp, utc, gl, gb)
    
    

if __name__ == '__main__':
    _main()
