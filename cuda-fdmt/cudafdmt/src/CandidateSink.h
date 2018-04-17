/*
 * CandidateSink.h
 *
 *  Created on: 28 Oct 2016
 *      Author: ban115
 */

#ifndef CANDIDATESINK_H_
#define CANDIDATESINK_H_

#include "DataSource.h"
#include <stdio.h>
#include <sys/socket.h>
#include <netinet/in.h>

class CandidateSink {
public:
	const size_t MAX_PACKET_BYTE = 1500;
	CandidateSink(DataSource* srcfile, const char* filename, bool negdm, const char* sockaddr, short sockprt);
	virtual ~CandidateSink();
	void add_candidate(int beam, int idt, int t, int ibc, float sn);
	void flush();// run after adding all candidates in a block


private:
	void flush_packet(); // flushes a packet
	DataSource* m_srcfile;
	FILE* m_candfile;
	bool m_negdm;
	int m_sockfd;
	struct sockaddr_in m_sockaddr;
	char* m_sockbuf;
	size_t m_sockbuf_size;
};

#endif /* CANDIDATESINK_H_ */
