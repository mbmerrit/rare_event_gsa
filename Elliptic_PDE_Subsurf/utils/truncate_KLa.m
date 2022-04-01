function KLa = truncate_KLa(KLa, ntrunc)

KLa.nKL = ntrunc;
KLa.sv = KLa.sv(1:ntrunc);
KLa.basis = KLa.basis(:,1:ntrunc);
