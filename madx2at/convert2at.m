clear
%% FCC lattice
fcc = 'z';
if fcc=='h'
    fname = '../fcc_lattice/Optics_FCCee_h_217_nosol_3.seq';
    seqfilemadX='../fcc_lattice/Optics_FCCee_h_217_nosol_3_for_at.seq';
    seqname = 'L000013';
    E0=120e9;
end
if fcc=='t'
    fname = '../fcc_lattice/Optics_FCCee_t_217_nosol_3.seq';
    seqfilemadX='../fcc_lattice/Optics_FCCee_t_217_nosol_3_for_at.seq';
    seqname = 'L000013';
    E0=175e9;
end
if fcc=='w'
    fname = '../fcc_lattice/Optics_FCCee_w_217_nosol_1.seq';
    seqfilemadX='../fcc_lattice/Optics_FCCee_w_217_nosol_1_for_at.seq';
    seqname = 'L000013';
    E0=80e9;
end
if fcc=='z'
    fname = '../fcc_lattice/Optics_FCCee_z_217_nosol_20.seq';
    seqfilemadX='../fcc_lattice/Optics_FCCee_z_217_nosol_20_for_at.seq';
    seqname = 'L000013';
    E0=45.6e9;
end
%t: ttbar, Ebeam = 175 GeV or 182.5 GeV h: ZH, Ebeam = 120 GeV w: W+W-, Ebeam = 80 GeV z: Z, Ebeam = 45.6 GeV
% execute conversion. file ..._AT_LATTICE.mat will be created.
mod = py.importlib.import_module('prepare_fcc_seq');
mod.reformat_seq(fname);
mod.save_seq(fname);
atfrommadx(seqfilemadX,E0);
!mv Optics* ../fcc_lattice/.
load([seqfilemadX(1:end-4) '_AT_LATTICE.mat']); 
eval(['ring = ',seqname,';']);
eval(['save ../fcc_lattice/fcc',fcc,'.mat ring']);
gen_lat_fun(strcat('../fcc_lattice/fcc',fcc))
!rm ../fcc_lattice/*for_at*