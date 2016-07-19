            if hdr.packetType != 129:
                continue

            #print hdr, address, len(data)#, binascii.hexlify(data)
            nint = int_cycles # get configured value, don't assume the packet size is right
            hdrsize = nint*4*8
            expected_size = hdrsize + 4*NBEAMS*NCHANS*nint
            actual_size = len(data)
            num_padbytes = expected_size - actual_size
            if expected_size != actual_size:
                print 'Expected size', expected_size, 'actual size', actual_size
            assert num_padbytes >= 0

            print hdr, 'BAT=', hex(hdr.bat), len(data), 'Is all zeros?', np.all(data==0), binascii.hexlify(data)

            if num_padbytes > 0:
                data += '\x00'*num_padbytes

            times = np.fromstring(data[0:hdrsize], dtype=np.uint64)
            lt = len(times)
            assert lt % 2 == 0
            frames = times[0:lt/2]
            bats = times[lt/2:]
            startFrames = frames[0::2]
            stopFrames = frames[1::2]
            startBats = bats[0::2]
            stopBats = bats[1::2]


            print '\t startBATs', ' '.join(map(hex, startBats))
            print '\t stopBATs', ' '.join(map(hex, stopBats))
            print '\t diffs', (stopBats - startBats)
            print '\t startFrames', startFrames
            print '\t stopFrames', stopFrames
            print '\t diffs', (stopFrames - startFrames)

            spectra = np.fromstring(data[hdrsize:], dtype=np.float32)
            # According to beamformer.cc the it's in FTB ordering, with
            # beam going fastest
            spectra.shape = (NCHANS, nint, NBEAMS)
            
            #plots.update(spectra)
            
#            n = (1184*2/8)*8
#            du = np.fromstring(data[0:n], dtype=np.uint64)
#            for i in xrange(len(du)):
#                print i, hex(du[i])


        #if spfile is None:
        #spfile = SigprocFile(values.output, 'w', hdr)

            #write_result_pkl(pkl_file, result)
        except:
            logging.exception('Exception processing spectra')
            break
            
    pkl_file.close()
    if bf is not None:
        del bf
