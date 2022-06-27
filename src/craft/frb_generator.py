

def gen_frbs():

    for i in range(nfrbs):
        blk = np.zeros(asdfasdfasdfads)

        yield blk
    


def gen_blocks():

    frb_gen = gen_frbs()

    frb_blk = next(frb_gen)

    while True:
        blk = np.random.rnandn(nbl, nc, nt)
        blk += frb_blk[: :, time_range]

        yield blk

        if finised_frb:
            frb_blk = next(frb_gen)
