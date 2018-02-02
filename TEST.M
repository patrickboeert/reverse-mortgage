function test

hecm = model();
sim  = simulation(hecm);
sim  = sim.get_demand([0.5,0.5]);

sim = sim.originate;
sim = sim.simulate;
sim = sim.close;
sim = sim.aggregate;

end