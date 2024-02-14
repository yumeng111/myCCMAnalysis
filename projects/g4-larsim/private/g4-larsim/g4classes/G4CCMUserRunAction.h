/*
   header file for G4CCMUserRunAction

   does not contain any actions.
*/

#ifndef g4_larsim_G4CCMUserRunAction_H
#define g4_larsim_G4CCMUserRunAction_H

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;

class G4CCMUserRunAction : public G4UserRunAction {
    public:
        G4CCMUserRunAction();
        virtual ~G4CCMUserRunAction();

        virtual void BeginOfRunAction(const G4Run* run);
        virtual void   EndOfRunAction(const G4Run* run);

    private:
};

#endif // g4_larsim_G4CCMUserRunAction_H
