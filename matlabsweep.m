% COMSOL / Matlab livelink script for minDE dispersed membrane model
% Physics @ Brandeis University 2023-2025, Michael M. Norton 
% 

% in order to save graphic outputs, use the following command when starting the comsol server: 
% $comsol server -graphics


datasetnamestr = 'sweep_name';

nD_range = [50, 200, 400, 800, 1000, 2000, 3000];
nE_range = [50, 200, 400, 800, 1000, 2000, 3000, 4000, 6000, 8000];
sim_time = 4000;

alpha_range = 1;
beta_range = 1;%10.^(linspace(-3,0,10));
gamma_range = 1; %10.^linspace(0,2,10);
paramcombs = combvec(nD_range,nE_range, alpha_range, beta_range, gamma_range);

[~, N_sims] = size(paramcombs);

for i_sweep = [1:16, 26:N_sims]%17:N_sims
    disp(['preparing simulation ' num2str(i_sweep) ' of ' num2str(N_sims)])
    nD= paramcombs(1,i_sweep);
    nE= paramcombs(2,i_sweep);
    alpha= paramcombs(3,i_sweep);
    beta= paramcombs(4,i_sweep);
    gamma = paramcombs(5,i_sweep);        

    disp(['  nD = ' num2str(nD)])
    disp(['  nE = ' num2str(nE)])
    disp(['  alpha = ' num2str(alpha)])
    disp(['  beta = ' num2str(beta)])
    disp(['  gamma = ' num2str(gamma)])


    filenamestr = ['sweep_',num2str(i_sweep)];
    pathstr = ['/savepath/', datasetnamestr,'/'];
    fullfilenamestr = [pathstr,filenamestr,'.mph'];
    fullfilenamestr_empty = [pathstr,filenamestr,'_empty.mph'];

    minEstr = [pathstr,'sweep_',num2str(i_sweep),'_images/minE_frames'];
    minDstr = [pathstr,'sweep_',num2str(i_sweep),'_images/minD_frames'];


    mkdir([pathstr,'sweep_',num2str(i_sweep),'_images'])
    mkdir(minEstr)
    mkdir(minDstr)

    import com.comsol.model.*
    import com.comsol.model.util.*

    model = ModelUtil.create('Model');

    model.modelPath('/modelpath');

    model.label('weakform_20240904_spiral.mph');

    model.comments(['Untitled\n\n']);

    model.param.set('alpha', num2str(alpha)); %1
    model.param.set('beta', num2str(beta)); %0.1
    model.param.set('gamma', num2str(gamma));
    model.param.set('kde', '0.34');
    model.param.set('ked', '0.01');
    model.param.set('kdD', '0.02');
    model.param.set('kD', '0.065');
    model.param.set('kdEr', '0.126');
    model.param.set('kdEl', '0.0002/10');
    model.param.set('ke', '0.0001*100*10');
    model.param.set('mu', '100'); %all minE remains active
    model.param.set('lamDD', '6');
    model.param.set('nD', num2str(nD));
    model.param.set('nE', num2str(nE));
    model.param.set('DiffuDD', 'DiffC');
    model.param.set('DiffuDT', 'DiffC');
    model.param.set('DiffuEr', 'DiffC');
    model.param.set('DiffuEl', 'DiffC');
    model.param.set('Diffud', '0.013*gamma');
    model.param.set('Diffude', '0.013*gamma');
    model.param.set('Diffue', '0.005*gamma');
    model.param.set('LX', '50');
    model.param.set('LY', '50');
    model.param.set('noisescale', '0.25');
    model.param.set('IC_amplitude', '0.0');
    model.param.set('DiffC', '60');
    model.param.set('Diffud_biharm', 'Diffud*lengthcutoff^2');
    model.param.set('Diffude_biharm', 'Diffude*lengthcutoff^2');
    model.param.set('Diffue_biharm', 'Diffue*lengthcutoff^2');
    model.param.set('lengthcutoff', '0.1');
    model.param.set('cmax_surf', num2str(2e4));
    model.param.set('cmax_vol', 'cmax_surf*alpha');

    model.modelNode.create('comp1');

    model.geom.create('geom1', 2);

    model.func.create('rn1', 'Random');
    model.func('rn1').set('type', 'normal');
    model.func('rn1').set('nargs', '2');

    model.mesh.create('mesh1', 'geom1');

    model.geom('geom1').create('r1', 'Rectangle');
    model.geom('geom1').feature('r1').set('size', {'LX' 'LY'});
    model.geom('geom1').feature('r1').set('base', 'center');
    model.geom('geom1').run;

    model.view.create('view2', 3);

    model.physics.create('w', 'WeakFormPDE', 'geom1');
    model.physics('w').field('dimensionless').component({'uDD' 'uDT' 'uEr' 'uEl' 'ud' 'ude' 'ue'});
    model.physics('w').create('wfeq2', 'WeakFormPDE', 2);
    model.physics('w').feature('wfeq2').selection.all;
    model.physics('w').create('wfeq3', 'WeakFormPDE', 2);
    model.physics('w').feature('wfeq3').selection.all;


    model.mesh('mesh1').create('auto_f1', 'FreeTri');



    model.physics('w').feature('wfeq1').set('weak', {'(-test(uDDx)*uDDx-test(uDDy)*uDDy)*DiffuDD'; '(-test(uDTx)*uDTx-test(uDTy)*uDTy)*DiffuDT'; '(-test(uErx)*uErx-test(uEry)*uEry)*DiffuEr'; '(-test(uElx)*uElx-test(uEly)*uEly)*DiffuEl'; '(-test(udx)*udx-test(udy)*udy)*Diffud'; '(-test(udex)*udex-test(udey)*udey)*Diffude'; '(-test(uex)*uex-test(uey)*uey)*Diffue'});
    model.physics('w').feature('wfeq1').label('diffusion');
    model.physics('w').feature('init1').set('uDD', '(nD/2)*(1+rn1(x,y)*noisescale+IC_amplitude*cos(2*pi*x/LX)*cos(2*pi*y/LY) )');
    model.physics('w').feature('init1').set('uEr', '(nE/2)*(1+rn1(x,y)*noisescale+IC_amplitude*cos(2*pi*x/LX)*cos(2*pi*y/LY) )');
    model.physics('w').feature('init1').set('ud', '(nD/2)*(1+rn1(x,y)*noisescale+IC_amplitude*cos(2*pi*x/LX)*cos(2*pi*y/LY) )');
    model.physics('w').feature('init1').set('ue', '(nE/2)*(1+rn1(x,y)*noisescale-IC_amplitude*cos(2*pi*x/LX)*cos(2*pi*y/LY))');

    model.physics('w').feature('wfeq2').set('weak',...
        {'(-lamDD*uDD+kde*ude)*test(uDD)';...
        '(lamDD*uDD-(kD*alpha*beta + beta*kdD*ud)*uDT*(cmax_vol-ud-ude)/cmax_vol)*test(uDT)';...
        '(-mu*uEr - ud*(beta*kdEr*uEr) + ke*ue)*test(uEr)';...
        '(mu*uEr - ud*(beta*kdEl*uEl))*test(uEl)';...
        '((beta*kD*alpha+ beta*kdD*ud)*uDT*(cmax_vol-ud-ude)/cmax_vol - ud*beta*(kdEr*uEr + kdEl*uEl) - ked*ud*ue/alpha)*test(ud)';...
        '(ud*beta*(kdEr*uEr + kdEl*uEl) - kde*ude   + ked*ud*ue/alpha)*test(ude)';...
        '(kde*ude - ked*ud*ue/alpha - ke*ue)*test(ue)'});
    model.physics('w').feature('wfeq2').label('reactions');
    model.physics('w').feature('wfeq3').set('weak', {'-uDDt*test(uDD)'; '-uDTt*test(uDT)'; '-uErt*test(uEr)'; '-uElt*test(uEl)'; '-udt*test(ud)'; '-udet*test(ude)'; '-uet*test(ue)'});
    model.physics('w').feature('wfeq3').label('time derivative');
    %model.physics('w').feature('wfeq6').set('weak', {'(-lamDD*uDD+kde*ude)*test(uDD)'; '(lamDD*uDD-(kD + kdD*ud)*uDT)*test(uDT)'; '(-mu*uEr - ud*(kdEr*uEr) + ke*ue)*test(uEr)'; '(mu*uEr - ud*(kdEl*uEl))*test(uEl)'; '((kD + kdD*ud)*uDT - ud*(kdEr*uEr + kdEl*uEl) - ked*ud*ue)*test(ud)'; '(ud*(kdEr*uEr + kdEl*uEl) - kde*ude   + ked*ud*ue)*test(ude)'; '(kde*ude - ked*ud*ue - ke*ue)*test(ue)'; '0'; '0'; '0'});
    %model.physics('w').feature('wfeq6').active(false);
    %model.physics('w').feature('wfeq6').label('reactions 1');

    model.mesh('mesh1').feature('size').set('custom', 'on');
    model.mesh('mesh1').feature('size').set('hmax', '1');
    model.mesh('mesh1').feature('size').set('hmin', '0.03');
    model.mesh('mesh1').run;

    % model.result.table('tbl1').comments('Surface Integration 1 ((uDD+uDT+ud+ude)/(LX*LY))');

    model.study.create('std1');
    model.study('std1').create('time', 'Transient');

    model.sol.create('sol1');
    model.sol('sol1').study('std1');
    model.sol('sol1').attach('std1');
    model.sol('sol1').create('st1', 'StudyStep');
    model.sol('sol1').create('v1', 'Variables');
    model.sol('sol1').create('t1', 'Time');
    model.sol('sol1').feature('t1').create('fc1', 'FullyCoupled');
    model.sol('sol1').feature('t1').create('st1', 'StopCondition');
    model.sol('sol1').feature('t1').feature.remove('fcDef');

    model.study('std1').feature('time').set('rtol', '1e-5');
    model.study('std1').feature('time').set('tlist', ['range(0,10,' num2str(sim_time) ')']);
    model.study('std1').feature('time').set('rtolactive', true);

    model.sol('sol1').attach('std1');
    model.sol('sol1').feature('t1').set('rtol', '1e-5');
    model.sol('sol1').feature('t1').feature('st1').set('stopconddesc', {'Stop expression 1'});
    model.sol('sol1').feature('t1').feature('st1').set('stopcondActive', {'on'});
    model.sol('sol1').feature('t1').feature('st1').set('stopcondterminateon', {'true'});
    model.sol('sol1').feature('t1').feature('st1').set('stopcondarr', {'timestep<1e-5'});
    model.sol('sol1').feature('t1').set('tlist', ['range(0,10,' num2str(sim_time) ')']);



    disp('saving empty model...')
    mphsave(model, fullfilenamestr_empty)
    disp('empty model saved!')


    
    try
        disp('running model...')
        model.sol('sol1').runAll;
        disp('run complete!')
        disp('saving model...')
        mphsave(model, fullfilenamestr)
        disp('model saved!')

        disp('generating images and movies...')

        model.result.create('pg1', 'PlotGroup2D');
        model.result.create('pg2', 'PlotGroup2D');
        model.result.create('pg3', 'PlotGroup2D');
        model.result('pg1').create('surf1', 'Surface');
        model.result('pg2').create('surf1', 'Surface');
        model.result('pg3').create('surf1', 'Surface');
        model.result.export.create('anim1', 'Animation');
        model.result.export.create('anim2', 'Animation');
        model.result.export.create('anim3', 'Animation');

        model.result('pg1').label('total minD');
        model.result('pg1').feature('surf1').set('descr', 'uDD+uDT+ud+ude');
        model.result('pg1').feature('surf1').set('colortable', 'GrayScale');
        model.result('pg1').feature('surf1').set('expr', 'uDD+uDT+ud+ude');
        model.result('pg2').label('total minE');
        model.result('pg2').feature('surf1').set('descr', 'uEr+uEl+ude+ue');
        model.result('pg2').feature('surf1').set('colortable', 'GrayScale');
        model.result('pg2').feature('surf1').set('expr', 'uEr+uEl+ude+ue');
        model.result('pg3').label('total minD green');
        model.result('pg3').feature('surf1').set('descr', 'uDD+uDT+ud+ude');
        model.result('pg3').feature('surf1').set('colortable', 'Green');
        model.result('pg3').feature('surf1').set('expr', 'uDD+uDT+ud+ude');
        model.result.export('anim1').label('exportframes minD');
        model.result.export('anim1').set('height', '800');
        model.result.export('anim1').set('framesel', 'all');
        model.result.export('anim1').set('width', '800');
        model.result.export('anim1').set('options', 'on');
        model.result.export('anim1').set('imagefilename', [minEstr,'/frames_.png']);
        model.result.export('anim1').set('type', 'imageseq');
        model.result.export('anim1').set('alwaysask', true);
        model.result.export('anim1').set('title', 'off');
        model.result.export('anim1').set('legend', 'off');
        model.result.export('anim1').set('logo', 'off');
        model.result.export('anim1').set('options', 'on');
        model.result.export('anim1').set('fontsize', '9');
        model.result.export('anim1').set('customcolor', [1 1 1]);
        model.result.export('anim1').set('background', 'color');
        model.result.export('anim1').set('axisorientation', 'on');
        model.result.export('anim1').set('grid', 'on');
        model.result.export('anim1').set('axes', 'off');
        model.result.export('anim2').label('exportframes minE');
        model.result.export('anim2').set('height', '800');
        model.result.export('anim2').set('framesel', 'all');
        model.result.export('anim2').set('width', '800');
        model.result.export('anim2').set('options', 'on');
        model.result.export('anim2').set('plotgroup', 'pg2');
        model.result.export('anim2').set('imagefilename', [minDstr,'/frames_.png']);
        model.result.export('anim2').set('type', 'imageseq');
        model.result.export('anim2').set('alwaysask', true);
        model.result.export('anim2').set('title', 'off');
        model.result.export('anim2').set('legend', 'off');
        model.result.export('anim2').set('logo', 'off');
        model.result.export('anim2').set('options', 'on');
        model.result.export('anim2').set('fontsize', '9');
        model.result.export('anim2').set('customcolor', [1 1 1]);
        model.result.export('anim2').set('background', 'color');
        model.result.export('anim2').set('axisorientation', 'on');
        model.result.export('anim2').set('grid', 'on');
        model.result.export('anim2').set('axes', 'off');




        model.result.export('anim3').label('exportgif minD green');
        model.result.export('anim3').set('giffilename', [pathstr,'sweep_',num2str(i_sweep),'_images/minD.gif']);
        model.result.export('anim3').set('plotgroup', 'pg3');
        model.result.export('anim3').set('title', 'off');
        model.result.export('anim3').set('legend', 'off');
        model.result.export('anim3').set('logo', 'off');
        model.result.export('anim3').set('options', 'off');
        model.result.export('anim3').set('fontsize', '9');
        model.result.export('anim3').set('customcolor', [1 1 1]);
        model.result.export('anim3').set('background', 'color');
        model.result.export('anim3').set('axisorientation', 'on');
        model.result.export('anim3').set('grid', 'on');
        model.result.export('anim3').set('axes', 'off');

        model.result.export('anim3').set('solnumtype', 'level1');
        model.result.export('anim3').set('framesel', 'all');
        model.result.export('anim1').run
        model.result.export('anim2').run
        model.result.export('anim3').run

        disp('done generating images and moives!')

    catch
        disp('model failed')
    end
        

    

    %out = model;


end
