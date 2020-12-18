function Rout = bb_kick(elemData,Rin)
kick_fun = elemData.lib.kick_calc;
spart = size(Rin);
Rout = Rin;
if spart>1
    for i = 1:spart(1)
        kick = kick_fun(Rin(i,1),Rin(i,3));
        Rout(i,2)=Rin(i,2)+kick{1};
        Rout(i,4)=Rin(i,4)+kick{2};
    end
else
    kick = kick_fun(Rin(1),Rin(3));
    Rout(2)=Rin(2)+kick{1};
    Rout(4)=Rin(4)+kick{2};
end

Rout


%amp = elemData.Amplitude;
%noise = ones(size(Rin)).*amp;
%Rout = Rin + noise;
end