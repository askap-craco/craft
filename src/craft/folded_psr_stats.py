#!/usr/bin/env python

from numpy import sqrt
import os
import logging
import glob
import subprocess as sb

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

    args = parser.parse_args()
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    influxout = None
    client = None
    if args.outfile:
        influxout = open(args.outfile, 'w')
    else:
        from influxdb import InfluxDBClient
        client = InfluxDBClient(host='akingest01', database='craft',
                                username='craftwriter', password='craft')
    body = []
    for filename in args.files:
        try:
            logging.debug("Handling {}".format(filename))
            body.append(get_folded_stats(filename, influxout,
                                         args, args.verbose))
        except Exception as inst:
            logging.exception('Exception in get_folded_stats:{0}'.format(inst))
    if influxout is not None:
        logging.debug("Writing data to file %s", str(args.outfile))
        for entry in body:
            influxout.write(str(entry)+"\n")
        logging.debug("wrote %s", str(body))
    if client is not None:
        logging.debug("Writing data to client %s", str(body))
        client.write_points(body)
    return 0


def get_basic_meta(filename_p, verbose):
    filename = os.path.basename(filename_p)
    filebits = filename.split('_')
    sbid = filebits[0]

    name_extra = ""
    if filename.endswith('_add.ar'):
        name_extra = "_ics"
        beam = filebits[2]
        coarse_timestamp = filebits[1]
        ant_count = len(glob.glob(os.path.dirname(filename_p)+"/"+sbid+'_' +
                                  coarse_timestamp+'_*_'+beam+'.ar'))
        ant = "ics_"+str(ant_count)
    else:
        ant = filebits[2]
        beam = filebits[3].split('.')[0]

    beam = beam[1:] if beam.startswith('0') else beam

    return sbid, name_extra, beam, ant


def get_folded_stats(filename_p, influxout, values, verbose):
    assert filename_p.endswith('.ar')
    sbid, name_extra, beam, ant = get_basic_meta(filename_p, verbose)

    name, tstamp, weff, on_count, off_rms, on_max,\
        off_avg, tint, snr, snr_pdmp, nchan, nsubint,\
        nchan_zapped = extract_snrs_from_archive(filename_p, verbose)
    snr_max = -1.0
    if float(off_rms) != 0.0:
        snr_max = (float(on_max)-float(off_avg)) / float(off_rms)
        snr_max = snr_max / sqrt(tint)

    snr = snr / sqrt(tint)
    snr_pdmp = snr_pdmp / sqrt(tint)

    zap_frac = float(nchan_zapped)/float(nchan*nsubint)
    snr_pdmp_zap = snr_pdmp / sqrt(1.-zap_frac)

    toa_unc, toa_gof, freq = get_timing_stats(filename_p, name, verbose)
    toa_unc = toa_unc * sqrt(tint)

    flux, flux_err, snr_flux = get_psr_flux(filename_p, freq, name, verbose)
    snr_flux = snr_flux / sqrt(tint)

    fields_dict = {'weff': weff, 'oncount': on_count,
                   'offrms': off_rms, 'onmax': on_max,
                   'offavg': off_avg, 'snr': snr,
                   'snr_max': snr_max, 'snr_pdmp': snr_pdmp,
                   'snr_pdmp_zap': snr_pdmp_zap, 'zap_frac': zap_frac,
                   'tint': tint, 'toa_unc': toa_unc, 'toa_gof': toa_gof,
                   'flux': flux, 'flux_err': flux_err, 'snr_flux': snr_flux}

    body = {'measurement': 'psrfold',
            'tags': {'psr': name+name_extra, 'sbid': sbid, 'ant': ant,
                     'beam': int(beam)},
            'time': int(tstamp*1e9),
            'fields': fields_dict}

    return body


def get_timing_stats(filename_p, name, verbose):
    pat_command = '/home/sha355/bin/pat -j FTp -a "/home/craftop/psr_templates/' + name\
                  + '*std" -A FDM -f "tempo2 IPTA" ' + filename_p
    if verbose:
        logging.debug("Running {}".format(pat_command))
    pat_process = sb.Popen(pat_command, shell=True, stdout=sb.PIPE,
                           stderr=sb.PIPE)
    (pat_out, pat_err) = pat_process.communicate()

    if pat_err:
        if "WARNING" in pat_err:
            logging.warn("WARNING while running pat: {}".format(pat_err))
        else:
            logging.error("Error while running pat: {}".format(pat_err))
            return -1.0, -1.0, ""
    pat_out_lines = pat_out.split('\n')
    if len(pat_out_lines) > 3:
        print("Error, expected 3 lines from pat, got {}: {}".format(
                len(pat_out_lines), pat_out))
    if "FORMAT" not in pat_out_lines[0]:
        print("Problem with pat output:", pat_out)
        return -1.0, -1.0, ""
    # This should:
    # file freq SAT unc tel fe_flag fe be_flag be f_flag f bw_flag bw
    # tobs_flag tobs tmplt_flag tmplt gof_flag gof nbin_flag nbin
    toa = pat_out_lines[1].split()
    toa_unc = float(toa[3])
    toa_gof = float(toa[18])
    freq_label = toa[16].split('_')[-1].split('.')[0]

    return toa_unc, toa_gof, freq_label


def get_psr_flux(filename_p, freq, name, verbose):
    psrflux_cmd = ['psrflux', '-j', 'FTp', '-a', '-s',
                   '/home/craftop/psr_templates/'+name+'_'+freq+'.std',
                   filename_p]
    if verbose:
        logging.debug("Running {}".format(psrflux_cmd))
    psrflux_proc = sb.Popen(psrflux_cmd, shell=False, stdout=sb.PIPE,
                            stderr=sb.PIPE)
    (psrflux_out, psrflux_err) = psrflux_proc.communicate()
    if psrflux_err:
        if "unloading" not in psrflux_err:
            print("Problem with psrflux:", psrflux_err)
            return -1.0, -1.0, -1.0
        psrflux_out_fn = filename_p + ".ds"
    with open(psrflux_out_fn) as fh:
        while True:
            line = fh.readline()
            line_el = line.split()
            if line_el[0] != '#':
                try:
                    flux = float(line_el[4])
                    flux_err = float(line_el[5])
                    snr_flux = flux / flux_err
                    return flux, flux_err, snr_flux
                except ValueError:
                    logging.error("Couldn't parse psrflux output for {}".format(filename_p))
                    logging.error(line)
                    return -1.0, -1.0, -1.0
                except ZeroDivisionError:
                    logging.error("Couldn't estimate S/N from psrflux for {}".format(filename_p))
                    logging.error("{} {}".format(flux, flux_err))
                    return flux, flux_err, -1.0
                except:
                    logging.error("Unknown error for {}".format(filename_p))
                    logging.error(line)
                    return -1.0, -1.0, -1.0
            if not line:
                break
    return -1.0, -1.0, -1.0


def extract_snrs_from_archive(archive_fn, verbose):
    from astropy.time import Time
    psrstat_command = ["/home/sha355/bin/psrstat", "-j", "DFTp", "-Qq", "-c",
                       "^off=minimum:smooth=mean:width=.5", "-c",
                       "name,int:mjd,weff,on:count,off:rms,on:max,off:avg,length,snr,snr=pdmp,snr",
                       archive_fn]
    if verbose:
        logging.debug("Running {}".format(psrstat_command))
    psrstat_process = sb.Popen(psrstat_command, shell=False, stdout=sb.PIPE,
                               stderr=sb.PIPE)
    (psrstat_out, psrstat_err) = psrstat_process.communicate()
    if psrstat_err:
        print("Error while running psrstat:", psrstat_err)
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
            tint = -1
            off_avg = -1
            print("Error interpreting psrstat output:")
            print(psrstat_out)
            print(len(out))

    weff_tmp = weff_turns.replace('.', '', 1).replace('e', '', 1).replace('-', '', 1)
    if not weff_tmp.isdigit():
        weff_turns = -1

    tstamp = Time(float(mjd), format="mjd").unix

    # get number of zapped channels:
    psrstat_command = ["/home/sha355/bin/psrstat", "-Qq", "-c",
                       "nchan,nsubint", archive_fn]
    if verbose:
        logging.debug("Running {}".format(psrstat_command))
    psrstat_process = sb.Popen(psrstat_command, shell=False, stdout=sb.PIPE,
                               stderr=sb.PIPE)
    (psrstat_out, psrstat_err) = psrstat_process.communicate()
    if verbose:
        logging.debug("Obtained {}".format(psrstat_out))
    out = psrstat_out.split()
    nchan = out[0]
    nsubint = out[1]

    psrstat_command = ["/home/sha355/bin/psrstat", "-Qq", "-c",
                       "int:wt", archive_fn]
    if verbose:
        logging.debug("Running {}".format(psrstat_command))
    psrstat_process = sb.Popen(psrstat_command, shell=False, stdout=sb.PIPE,
                               stderr=sb.PIPE)
    (psrstat_out, psrstat_err) = psrstat_process.communicate()
    out = psrstat_out.split()[0].split(",")
    nchan_zapped = out.count('0')
    if verbose:
        logging.debug("Determined {} zapped channels (out of {})".format(nchan_zapped, int(nchan)*int(nsubint)))

    if verbose:
        debug_str = "Obtained {0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {11} {12}"
        debug_str = debug_str.format(name, mjd, weff_turns, on_count,
                                     off_rms, on_max, off_avg, tint, snr, snr_pdmp,
                                     nchan, nsubint, nchan_zapped)
        logging.debug(debug_str)

    return name, tstamp, float(weff_turns), int(on_count), float(off_rms),\
        float(on_max), float(off_avg), float(tint), float(snr),\
        float(snr_pdmp), int(nchan), int(nsubint), int(nchan_zapped)


if __name__ == "__main__":
    _main()
