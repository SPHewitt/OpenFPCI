// Hack: Cannot write to timestep before 
// runTime.write() is called in the top level solver
// Write is done on subsequent timestep

//    if (runTime().outputTime())
    if(writePFlux_)
    {
        OFstream gradPFile
        (
            runTime().path()/writeTime_/"uniform"/"gradPf.raw"
        );

        if(gradPFile.good())
        {
            gradPFile << gradPfOld_ << endl;
        }
        else
        {
            //FatalErrorIn(args.executable())
            Info << "Cannot open file "
                << runTime().path()/writeTime_/"uniform"/"gradP.raw"
                << exit(FatalError);
        };
	writePFlux_ = false;
    };


    if(runTime().outputTime())
    {
	writePFlux_ = true;
	writeTime_ = runTime().timeName();
	gradPfOld_ = gradPf_;
    }   
