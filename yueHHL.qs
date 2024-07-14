// HHL Algorithm
//
//#Description
//Used to solve linear systems of equations
//In this implementation:
// 1. **PrepareState**: Prepares the initial state based on the coefficients.
// 2. **ApplyHamiltonianEvolution**: Applies the unitary operator representing the Hamiltonian evolution.
// 3. **QuantumPhaseEstimation**: Performs quantum phase estimation to find the eigenvalues of the unitary.
// 4. **ControlledRotation**: Applies controlled rotations based on the eigenvalues.
// 5. **HHL**: Combines these components to implement the HHL algorithm

namespace HHL {
    open Microsoft.Quantum.Arrays;
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;
    open Microsoft.Quantum.Unstable.StatePreparation;

    @EntryPoint()
    operation Main() : Unit {

        let scale = 2;
        let vector = [1.0 / Sqrt(2.0), 1.0 / Sqrt(2.0)];
        let matrix = [
            [1.0, 0.0],
            [0.0, - 1.0]
        ];

        use stateVectorb = Qubit[1];
        PreparePureStateD(vector, stateVectorb);

        mutable result = false;
        repeat {
            set result = ApplyHHL(HamiltonianEvolutionSample, stateVectorb);
        } until result;
        
        DumpRegister(stateVectorb);
        Reset(stateVectorb[0]);
    }

    operation ApplyHHL(unitaryA : Qubit[] => Unit is Adj + Ctl, targetRegister : Qubit[]) : Bool {

        let numClockQubits = CalculateNumClockQubits();
        let evolutionTime = CalculateEvolutionTime();
        // let unitary = HamiltonianEvolution(A, evolutionTime, targetRegister);
        use clockRegister = Qubit[numClockQubits];
        use anciliaRegister = Qubit();

        within {
            // Unitary = HamiltonianSimulation?
            PhaseEstimation(unitaryA, clockRegister, targetRegister);
        } apply {
            let scaling = CalculateScaling();
            let negVal = true;
            Reciprocal(scaling, negVal, clockRegister, anciliaRegister);
        }

        ResetAll(clockRegister);
        let postSelectResult = M(anciliaRegister);
        Reset(anciliaRegister);

        if postSelectResult == One {
            return true;
        } else {
            return false;
        }
    }


    internal operation CalculateNumClockQubits() : Int {
        return 3; // implement later
    }

    internal operation CalculateEvolutionTime() : Double {
        return PI() / 2.0; // implement later
    }

    internal operation CalculateScaling() : Double {
        return 0.25; // implement later
    }


    operation PrepareUniform(inputQubits : Qubit[]) : Unit is Adj + Ctl {
        for q in inputQubits {
            H(q);
        }
    }

    operation PhaseEstimation(unitary : (Qubit[] => Unit is Adj + Ctl), clockQubits : Qubit[], phiQubits : Qubit[]) : Unit is Adj + Ctl {

        // Apply Hadamard gates to all clock qubits
        PrepareUniform(clockQubits);

        // Apply controlled unitary operations
        let nClock = Length(clockQubits);
        for i in 0..nClock - 1 {
            let power = 2^i; // little-endian, first qubits present less significant bits
            for _ in 0..power - 1 {
                Controlled unitary([clockQubits[i]], phiQubits);
            }
        }

        // Apply inverse Quantum Fourier Transform
        Adjoint ApplyQFT(clockQubits);
    }

    operation Reciprocal(scaling : Double, negVal : Bool, clockQubits : Qubit[], anciliaQubit : Qubit) : Unit {
        mutable negValInt = 0;
        if negVal {
            set negValInt = 1;
        }

        let nClock = Length(clockQubits);
        let nAbsClock = nClock - negValInt;

        let nVal = 2^nClock;
        let nAbsVal = 2^nAbsClock;

        for i in 1..nAbsVal- 1 {
            let angle = 2.0 * ArcSin(scaling * IntAsDouble(nAbsVal) / IntAsDouble(i));
            // Message($"arcsin(x), x = {scaling * IntAsDouble(nAbsVal) / IntAsDouble(i)};angle = {angle}; bitsize i = {i}");
            ApplyControlledOnInt(i, Ry(angle, _), clockQubits[0..nAbsClock - 1], anciliaQubit);
        }

        if negVal {
            // Message("here in negVal");
            for i in 1..nAbsVal -1 { // counteract
                let negAngle = - 2.0 * ArcSin(scaling * IntAsDouble(nAbsVal) / IntAsDouble(i));
                // Message($"arcsin(x), x = {scaling * IntAsDouble(nAbsVal) / IntAsDouble(i)};negAngle = {negAngle}; bitsize i = {i}");
                Controlled ApplyControlledOnInt([Tail(clockQubits)], (i, Ry(negAngle, _), clockQubits[0..nAbsClock - 1], anciliaQubit));
            }

            for i in 1..nAbsVal-1 { // two's complement representation
                let negAngle = 2.0 * ArcSin(scaling * (- 1.0) / (1.0 -IntAsDouble(i) / IntAsDouble(nAbsVal)) );
                // Message($"arcsin(x), x = {scaling * (- 1.0) / (1.0 -IntAsDouble(i) / IntAsDouble(nAbsVal))} ; negAngle = {negAngle}; bitsize i = {i}");
                Controlled ApplyControlledOnInt([Tail(clockQubits)], (i, Ry(negAngle, _), clockQubits[0..nAbsClock - 1], anciliaQubit));
            }

        }

    }

    operation HamiltonianEvolutionSample(targetRegister : Qubit[]) : Unit is Adj + Ctl {

        Y(targetRegister[0]);
        X(targetRegister[0]);
    }
}