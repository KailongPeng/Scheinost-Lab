function [mask,bigval] = lagestcc(a)
%Brett's shot at Doug's "Biggest 4-connected" puzzler
%08/18/08
lindices = tril(ones(268));
lindices(find(eye(268)))=0;
lindices = find(lindices);

network = zeros(268,268);
network(lindices) = a;
network = network'+network;

biggestArea = 0;
L = bwlabel(network==1,4);
stats = regionprops(L,'area');
[tmp,tmpbig] = max([stats.Area]);
if tmp > biggestArea
    biggestArea = tmp;
    mask = L == tmpbig;
    bigval = 1;
end
end