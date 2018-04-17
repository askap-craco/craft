/*
 * CandidateSink.cpp
 *
 *  Created on: 28 Oct 2016
 *      Author: ban115
 */

#include "CandidateSink.h"
#include <stdio.h>
#include <stdlib.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <assert.h>

CandidateSink::CandidateSink(DataSource* srcfile, const char* filename,
		bool negdm, const char* sockaddr, short sockprt) {
	m_srcfile = srcfile;
	m_negdm = negdm;
	if (filename == NULL) {
		char fname[4096];
		sprintf(fname, "%s.cand", srcfile->name());
		m_candfile = fopen(fname, "w");
	} else {
		m_candfile = fopen(filename, "w");
	}
	if (!m_candfile) {
		perror("Could not open candidate file");
		exit(EXIT_FAILURE);
	}
	fprintf(m_candfile, "# S/N, sampno, secs from file start, boxcar, idt, dm, beamno\n");

	if (sockaddr != NULL && strlen(sockaddr) > 0 && sockprt > 0) {
		if ((m_sockfd = socket(AF_INET, SOCK_DGRAM, 0)) < 0) {
			perror("Cannot create socket");
			exit(EXIT_FAILURE);
		}
		struct hostent* hostname = gethostbyname(sockaddr);
		if (hostname  == NULL) {
			herror("Invalid socket address");
			exit(EXIT_FAILURE);
		}
		struct in_addr **addr_list = (struct in_addr **) hostname->h_addr_list;
		assert(addr_list[0] != NULL);
		m_sockaddr.sin_addr = *addr_list[0];
		m_sockaddr.sin_family = AF_INET;
		m_sockaddr.sin_port = htons(sockprt);
		m_sockbuf = (char*) malloc(MAX_PACKET_BYTE);
		m_sockbuf_size = 0;
	} else {
		m_sockbuf = NULL;
		m_sockbuf_size = 0;
	}
}

CandidateSink::~CandidateSink() {
	if (m_candfile) {
		fflush(m_candfile);
		fclose(m_candfile);
	}
	if (m_sockbuf) {
		free(m_sockbuf);
		// TODO: release socket
	}
}

void CandidateSink::add_candidate(int ibeam, int idt, int t, int ibc, float sn)
{
	// Typical line: 11.1591 7505890 9499.64 12      178     736.883 6       7505720 7507042
	// S/N sample_number dunno width dunno dm dunno dunno dunno
	// dt = 4.15ms * DM * (nu1**-2 - nu2**-2)
	// DM =

	float  dm = m_srcfile->dm_of_idt(idt);
	double time_from_file = t*m_srcfile->tsamp();
	if (m_negdm) {
		idt = -idt;
		dm = -dm;
	}
	fprintf(m_candfile, "%f %lu %f %d %d %0.3f %d\n", sn, t, time_from_file,
			ibc, idt, dm, ibeam);

	if (m_sockbuf != NULL) {
		int line_len = snprintf(m_sockbuf+m_sockbuf_size, MAX_PACKET_BYTE, "%f %lu %f %d %d %0.3f %d\n", sn, t, time_from_file,
				ibc, idt, dm, ibeam);
		m_sockbuf_size += line_len;

		// flush packet early if it's getting close to being big
		if (m_sockbuf_size + 2*line_len >= MAX_PACKET_BYTE) {
			flush_packet();
		}
	}
}

void CandidateSink::flush_packet()
{
	if (m_sockbuf != NULL && m_sockbuf_size > 0) {
		if(sendto(m_sockfd, m_sockbuf, strlen(m_sockbuf), 0,
				(struct sockaddr*)&m_sockaddr,
				sizeof(struct sockaddr_in)) < 0) {
			printf("Could not send data to UDP socket\n");
		}
		m_sockbuf[0] = 0;
		m_sockbuf_size = 0;
	}
}

void CandidateSink::flush()
{
	fflush(m_candfile);
	flush_packet();
}

