function g=CVgroups(Ngroup,bin_number)

g=zeros(size(bin_number));
Nbin=max(bin_number);
for i=1:Nbin
    bin=bin_number==i;
    g(bin)=mod(randperm(sum(bin)),Ngroup);
end

g(g==0)=Ngroup;