zmq_setsockopt(socket, ZMQ_NOBLOCK);

if (zmq_recv(socket, message) == 0) {
        // Successful receive
        // Apply calibration
} else {
    if (errno == EAGAIN) {
            // there was no message to receive
            } else {
                // There was a big bad problem.
            }
    
}

struct _uvaver_header_t {
    uint64_t timestamp;
    size_t num_baselines;
    size_t num_pol;
    size_t integration_time ; // other meta data
    
    
} uvaver_header_t;

header followed by nbaliens * npol ddata values

or we cave protocol defines 


while  True:
    hdr_bytes = zmqsocket.read(ZMQ_RCVMORE)
    data_bytes = zmqsocket.read()
    hdr = protocol_class.load(data[0:procotol_class.size])
    data = np.frombytes(data[protocol_class.size:], dtype=np.complex64).reshape(hdr.num_baselines, hdr.num_pol)
