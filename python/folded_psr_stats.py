#!/usr/bin/env python

from numpy import sqrt
import os
import logging

__author__ = "Stefan Oslowski <stefanoslowski@swin.edu.au>"


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description',
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
                        help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.add_argument('-O', '--outfile',
                        help='Output file for influxdb data. Inhibits'
                        ' live loading')

    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    influxout = None
    client = None
    if values.outfile:
        influxout = open(values.outfile, 'w')
    else:
        from influxdb import InfluxDBClient
        client = InfluxDBClient(host='akingest01', database='craft')
        # print "not available yet yo"
        # sys.exit(0)
    body = []
    for filename in values.files:
        try:
            logging.debug("Handling {}".format(filename))
            body.append(get_folded_stats(filename, client, influxout, values))
        except Exception as inst:
            logging.exception('Blah exception in get_folded_stats: {0}'.format(inst))
    if influxout is not None:
        logging.debug("Writing data to file %s", str(values.outfile))
        for entry in body:
            influxout.write(str(entry)+"\n")
        logging.debug("wrote %s", str(body))
    if client is not None:
        logging.debug("Writing data to client %s", str(body))
        client.write_points(body)
    return 0


def get_folded_stats(filename_p, client, influxout, values):
    assert filename_p.endswith('.ar')
    filename = os.path.basename(filename_p)
    filebits = filename.split('_')
    sbid = filebits[0]
    # coarse_timestamp = filebits[1]
    ant = filebits[2]
    beam = filebits[3].split('.')[0]
    beam = beam[1:] if beam.startswith('0') else beam

    name, tstamp, weff, on_count, off_rms, on_max,\
        off_avg, tint, snr, snr_pdmp = extract_stats_from_archive(filename_p)
    snr_max = -1.0
    if float(off_rms) != 0.0:
        snr_max = (float(on_max)-float(off_avg)) / float(off_rms)
        snr_max = snr_max / sqrt(tint)

    snr = snr / sqrt(tint)
    snr_pdmp = snr_pdmp / sqrt(tint)

    fields_dict = {'weff': float(weff), 'oncount': int(on_count),
                   'offrms': float(off_rms), 'onmax': float(on_max),
                   'offavg': float(off_avg), 'snr': float(snr),
                   'snr_max': float(snr_max), 'snr_pdmp': float(snr_pdmp),
                   'tint': float(tint)}

    body = {'measurement': 'psrfold',
            'tags': {'psr': name, 'sbid': sbid, 'ant': ant, 'beam': int(beam)},
            'time': int(tstamp*1e9),
            'fields': fields_dict}

    return body


def extract_stats_from_archive(archive_fn):
    import subprocess as sb
    from astropy.time import Time
    psrstat_command = ["/home/sha355/bin/psrstat", "-j", "DFTp", "-Qq", "-c",
                       "name,int:mjd,weff,on:count,off:rms,on:max,off:avg,length,snr,snr=pdmp,snr",
                       archive_fn]
    psrstat_process = sb.Popen(psrstat_command, shell=False, stdout=sb.PIPE,
                               stderr=sb.PIPE)
    (psrstat_out, psrstat_err) = psrstat_process.communicate()
    if psrstat_err:
        print "Error while running psrstat:", psrstat_err
        return -1
    else:
        out = psrstat_out.split()
        if len(out) == 16:
            # most likely pdmp s/n failed
            name = out[0]
            mjd = out[1]
            weff_turns = out[2]
            on_count = out[3]
            off_rms = out[4]
            on_max = out[5]
            off_avg = out[6]
            tint = out[7]
            snr = out[8]
            snr_pdmp = -1
        elif len(out) == 10:
            name = out[0]
            mjd = out[1]
            weff_turns = out[2]
            on_count = out[3]
            off_rms = out[4]
            on_max = out[5]
            off_avg = out[6]
            tint = out[7]
            snr = out[8]
            snr_pdmp = out[9]
        else:
            # unknown error
            name = "Error"
            mjd = 55000
            weff_turns = -1
            on_count = -1
            off_rms = -1
            on_max = -1
            snr = -1
            snr_pdmp = -1
            off_avg = -1
            print "Error interpreting psrstat output:"
            print psrstat_out

    weff_tmp = weff_turns.replace('.', '', 1).replace('e', '', 1).replace('-', '', 1)

    if not weff_tmp.isdigit():
        weff_turns = -1

    tstamp = Time(float(mjd), format="mjd").unix
    return name, tstamp, float(weff_turns), float(on_count), float(off_rms),\
        float(on_max), float(off_avg), float(tint), float(snr),\
        float(snr_pdmp)


if __name__ == "__main__":
    _main()
