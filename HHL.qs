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
// 
namespace HHL {
    open Microsoft.Quantum.Diagnostics;
    open Microsoft.Quantum.Intrinsic;
    open Microsoft.Quantum.Canon;
    open Microsoft.Quantum.Math;
    open Microsoft.Quantum.Convert;

    @EntryPoint()
    operation EntryPointMain() : Unit {
        let applyControlledUnitary = ApplyControlledUnitary;
        let precisionQubits = 3;
        let coefficients = [1.0, 0.0, 0.0, 0.0];
        
        Main(applyControlledUnitary, precisionQubits, coefficients);
    }

    operation Main(ApplyControlledUnitary: (Qubit, Qubit[], Double) => Unit, precisionQubits: Int, coefficients: Double[]) : Unit {
        use targetRegister = Qubit[4];
        
        // Prepare |b⟩ state
        PrepareBState(targetRegister, coefficients);
        
        // Apply Quantum Phase Estimation with controlled rotation for eigenvalue inversion
        QuantumPhaseEstimationWithInversion(ApplyControlledUnitary, precisionQubits, targetRegister);
        
        // Display the state (for demonstration, requires measurement)
        DumpRegister(targetRegister);

        // Measure the final state of the target register
        for i in 0 .. Length(targetRegister) - 1 {
            let result = M(targetRegister[i]);
            Message($"Qubit {i}: {result}");
        }

        // Reset the qubits
        ResetAll(targetRegister);
    }

    function CalculateNorm(coefficients: Double[]): Double {
        mutable norm = 0.0;
        for i in 0 .. Length(coefficients) - 1 {
            set norm += coefficients[i] * coefficients[i];
        }
        return Sqrt(norm);
    }

    function NormalizeCoefficients(coefficients: Double[], norm: Double): Double[] {
        mutable normalizedCoefficients = coefficients;
        for i in 0 .. Length(coefficients) - 1 {
            set normalizedCoefficients w/= i <- coefficients[i] / norm;
        }
        return normalizedCoefficients;
    }

    operation PrepareBState(register: Qubit[], coefficients: Double[]): Unit is Adj + Ctl {
        let n = Length(register);
        if (n != Length(coefficients)) {
            fail "The length of coefficients array must match the number of qubits.";
        }

        // Normalize the coefficients using a function
        let norm = CalculateNorm(coefficients);
        let normalizedCoefficients = NormalizeCoefficients(coefficients, norm);

        // Prepare |b⟩ state using the normalized coefficients
        for i in 0 .. (n - 1) {
            Ry(2.0 * ArcCos(normalizedCoefficients[i]), register[i]);
        }
    }

    function PowI (a : Int, power : Int) : Int {
        mutable result = 1;
        for _ in 1..power {
            set result *= a;
        }
        return result;
    }

    function RoundAsInt(value : Double) : Int {
        return (Floor(value + 0.5));
    }

    operation HamiltonianSimulation(hamiltonian: Double[][], time: Double, qubits: Qubit[]): Unit is Adj + Ctl {
        let nQubits = Length(qubits);
        let nMatrixRows = Length(hamiltonian);
        let nMatrixCols = Length(hamiltonian[0]);
    
        // Basic validation
        if nMatrixRows != nMatrixCols {
            fail "Hamiltonian matrix must be square.";
        }
        if nMatrixRows != PowI(2, nQubits) {
            fail "Hamiltonian matrix dimensions do not match the number of qubits.";
        }

        // Time evolution parameters
        let stepSize = 0.1; // Time step for the Trotter-Suzuki decomposition
        let nSteps = RoundAsInt(time / stepSize);
        
        for _ in 1..nSteps {
            // Apply Hamiltonian simulation for each step
            for i in 0..nMatrixRows - 1 {
                for j in 0..nMatrixCols - 1 {
                    let h_ij = hamiltonian[i][j];
                    if (h_ij != 0.0) {
                        // Apply controlled Rz rotations based on the Hamiltonian matrix elements
                        // Adjust angles based on the time step size
                        let angle = -2.0 * h_ij * stepSize;
                        if (i == j) {
                            for k in 0..nQubits - 1 {
                                Rz(angle, qubits[k]);
                            }
                        }
                    }
                }
            }
        }
    }

    operation ApplyControlledUnitary(control: Qubit, register: Qubit[], time: Double): Unit is Adj + Ctl {
        // Example Hamiltonian matrix A (identity matrix for simplicity)
        let hamiltonian = [
            [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
        ];

        // Apply controlled Hamiltonian simulation
        let controlledHamiltonianSimulation = Controlled HamiltonianSimulation;
        controlledHamiltonianSimulation([control], (hamiltonian, time, register));
    }

    function PowD(x: Double, y: Double): Double {
        return (x^y);
    }

    operation QuantumPhaseEstimationWithInversion(
        ApplyControlledUnitary: (Qubit, Qubit[], Double) => Unit,
        precisionQubits: Int,
        targetRegister: Qubit[]
    ) : Unit {
        use precisionRegister = Qubit[precisionQubits];
        use ancilla = Qubit();
        
        within {
            for j in 0 .. precisionQubits - 1 {
                Controlled R1([precisionRegister[j]], (-IntAsDouble(j + 1), ancilla));
            }
            X(ancilla);
            H(ancilla);
            for j in 0 .. precisionQubits - 1 {
                Controlled R1([precisionRegister[j]], (IntAsDouble(j + 1), ancilla));
            }
            H(ancilla);
            X(ancilla);
        } apply {
            // Apply the rotation to the target qubits
            for j in 0 .. precisionQubits - 1 {
                let angle = 2.0 * PI() / PowD(2.0, IntAsDouble(j + 1));
                Controlled Rz([precisionRegister[j]], (angle, targetRegister[0]));
            }

            // Apply the provided controlled unitary operation
            for j in 0 .. precisionQubits - 1 {
                ApplyControlledUnitary(precisionRegister[j], targetRegister, IntAsDouble(j + 1));
            }
        }

        // Reset the precision and ancilla qubits
        ResetAll(precisionRegister);
        Reset(ancilla);
    }
}
