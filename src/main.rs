use chem_rs::{basic::*, core::*};

fn main() {
    let chain = Chain::new(Atom::Nitrogen)
        .add_atoms(Atom::Carbon, BondType::Single, 4)
        .convert_chain_type(ChainType::MonoCyclic)
        .add_side_chain(&oxygen(), 2, BondType::Double(false))
        .add_side_chain(&oxygen(), 5, BondType::Double(false));

    println!(
        "Chemfig: {}\nSMILES: {}",
        chain.chemfig_repr(),
        chain.smiles()
    );
}
