include ../config

ROOTDIR=..
LIBDIR=$(ROOTDIR)/lib
INCLUDE=-I$(ROOTDIR)/include -I$(PREFIX)/include

CC=gcc
CFLAGS=-ansi -Wall -fPIC -O3 $(INCLUDE) -funroll-loops
LD=gcc
LDFLAGS=-shared
SIGNUP=echo "ZM ver."$(VERSION)" Copyright (C) 1998 Tomomichi Sugihara (Zhidao)" >>

OBJ=zm_ieee.o zm_misc.o\
	zm_sf_erf.o zm_sf_gamma.o zm_sf_bessel.o\
	zm_rand.o zm_stat.o\
	zm_complex.o zm_complex_arith.o zm_complex_pe.o\
	zm_raw_vec.o zm_raw_mat.o\
	zm_vec.o zm_vec_array.o zm_vec_list.o zm_vec_tree.o zm_vec_ring.o zm_mat.o\
	zm_cvec.o zm_cmat.o\
	zm_le.o zm_le_pivot.o zm_le_lu.o zm_le_lq.o zm_le_minv.o zm_le_mpinv.o zm_le_tridiag.o zm_le_gen.o zm_le_lyapnov.o\
	zm_eig.o\
	zm_mca.o zm_mca_cluster.o zm_mca_gmm.o\
	zm_seq.o\
	zm_pex.o zm_pex_eq.o\
	zm_ip.o zm_ip_lin.o zm_ip_lagrange.o zm_ip_spline.o zm_ip_akima.o zm_ip_pex.o zm_ip_ipio.o\
	zm_nurbs.o\
	zm_opt_line.o\
	zm_opt_lp_stdcnv.o zm_opt_lp_simplex.o zm_opt_lp_pdip.o\
	zm_opt_lcp_lemke.o zm_opt_lcp_ip.o\
	zm_opt_qp.o zm_opt_nm.o zm_opt_ga.o zm_opt_dm.o\
	zm_nle_se.o zm_nle_dm.o zm_nle_ss.o\
	zm_ode.o zm_ode_dc.o zm_ode_euler.o zm_ode_heun.o zm_ode_rk4.o zm_ode_rkg.o zm_ode_erk.o zm_ode_adams.o zm_ode_beuler.o zm_ode_bk4.o zm_ode_gear.o zm_ode2.o\
	zm_fft.o\
	zm_intg.o zm_oscil_vdp.o zm_oscil_kura.o\
	zm_graph.o zm_graph_search.o\
	zm_rrt.o\
	zm_parse.o
DLIB=libzm.so

TARGET=$(DLIB)

all: $(TARGET)

$(DLIB): $(OBJ)
	@echo " LD	" $^
	@$(LD) $(LDFLAGS) -o $@ $^ > /dev/null; $(SIGNUP) $@
	-@mv $@ $(LIBDIR)
%.o: %.c
	@echo " CC	" $<
	@$(CC) $(CFLAGS) -c $< > /dev/null
clean:
	@echo " CLEAN"
	-@rm -f *.o *~ \#* $(TARGET)
