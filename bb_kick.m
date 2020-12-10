function Rout = bb_kick(elemData,Rin)
amp = elemData.Amplitude;
noise = ones(size(Rin)).*amp;
Rout = Rin + noise;
end