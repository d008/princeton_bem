%---- This is an example script for running a BEM simulation ----%

%% -- First: select a rotor from those available by creating a class object-- %%

rt = pbem_cls;
%User will be prompted to choose a rotor from those available

%% -- Second: set the run conditions -- %%
rt = rt.runsetup;
%NOTE: the run conditions can also be directly altered via dot notation 
%calls, i.e. to change the blade number to 4 enter: rt.runcon.nb = 4;

%% -- Third: Run the simulation! -- %%
rt = rt.run_pbem;
%Note: output is given in two structures, "bemd" which contains global
%results and "bld" which gives blade-level information

%% -- Fourth: Plot the results -- %%
hs = plot_g(bemd);
%NOTE: will give an overview plot of the global variables in bemd.