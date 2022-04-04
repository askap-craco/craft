
import fitswriter
import astropy.io.fits as fits
import numpy as np
import os

example_dtype = [('frame_id', '<u8'), # Words 1,2
                ('bat', '<u8'), # Words 3,4
                 ('beam_number','<u1'), # Word 5
                 ('sample_number','<u1'),
                 ('channel_number','<u1'),
                 ('fpga_id','<u1'), 
                 ('nprod','<u2'),# Word 6
                 ('flags','<u1'), # (sync_bit, last_packet, enable_debug_header, enable_test_data
                ('zero1','<u1'), # padding
                 ('zero2','<u4'), # word 7
                ('zero3','<u4'), # word 8
                ('data', '<i2',(2,2)) # data
]
   

def compare(npdata, fitsdata, byteswap=False):
    for k in npdata.dtype.fields.keys():
        d = np.array(fitsdata[k])
        if byteswap:
            d= d.byteswap()
        
        assert np.all(npdata[k] == d)

    
def test_write():
    w =  fitswriter.FitsTableWriter('test.fits', example_dtype)

    data = []
    for i in range(10):
        v = np.ones(1, dtype=example_dtype)
        data.append(v)
        w.write(v)

    w.close()

    f = fits.open('test.fits')
    print(f.info())
    assert len(f) == 2 # SHould have 2 HDUs

    for i in range(10):
        v = f[1].data[i]
        print(type(v), type(data[i]))
        compare(data[i],v)

    f.close()
    os.remove('test.fits')

def test_write_with():
    with fitswriter.FitsTableWriter('test.fits', example_dtype) as w:
        data = []
        for i in range(10):
            v = np.ones(1, dtype=example_dtype)
            data.append(v)
            print(v)
            w.write(v)
            print(v)

    f = fits.open('test.fits')
    print(f.info())
    assert len(f) == 2 # SHould have 2 HDUs

    for i in range(10):
        v = f[1].data[i]
        print(type(v), type(data[i]))
        compare(data[i],v)

    f.close()
    os.remove('test.fits')

def test_write_with_no_byteswap():
    with fitswriter.FitsTableWriter('test.fits', example_dtype, byteswap=False) as w:
        data = []
        for i in range(10):
            v = np.ones(1, dtype=example_dtype)
            data.append(v)
            print(v)
            w.write(v)
            print(v)

    f = fits.open('test.fits')
    print(f.info())
    assert len(f) == 2 # SHould have 2 HDUs

    for i in range(10):
        v = f[1].data[i]
        print(type(v), type(data[i]))
        compare(data[i],v, True)

    f.close()
    os.remove('test.fits')
    


            

    
    
                 
