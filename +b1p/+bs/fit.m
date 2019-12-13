function model = fit(dat,model,opt)

model = b1p.bs.mean.init(dat,model,opt);
model = b1p.bs.sens.init(dat,model,opt);
model = b1p.bs.b1.init(dat,model,opt);

for it=1:opt.nbiter
    
    model = b1p.bs.sens.update(dat,model,opt);
    model = b1p.bs.mean.update(dat,model,opt);
    model = b1p.bs.b1.update(dat,model,opt);
    
end