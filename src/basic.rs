//! This module contains some basic common building block compounds.

use crate::core::*;

/// Returns a carbon chain with length specified by the IUPAC root word.
/// Note: All bonds are single bonds and chain is open.
pub fn root_word_chain(word: &str) -> Chain {
    match word {
        "meth" => Chain::new(Atom::Carbon),
        "eth" => Chain::new_carbon_chain(2, ChainType::Open),
        "prop" => Chain::new_carbon_chain(3, ChainType::Open),
        "but" => Chain::new_carbon_chain(4, ChainType::Open),
        "pent" => Chain::new_carbon_chain(5, ChainType::Open),
        "hex" => Chain::new_carbon_chain(6, ChainType::Open),
        "hept" => Chain::new_carbon_chain(7, ChainType::Open),
        "oct" => Chain::new_carbon_chain(8, ChainType::Open),
        "non" => Chain::new_carbon_chain(9, ChainType::Open),
        "dec" => Chain::new_carbon_chain(10, ChainType::Open),
        "undec" => Chain::new_carbon_chain(11, ChainType::Open),
        "dodec" => Chain::new_carbon_chain(12, ChainType::Open),
        "tridec" => Chain::new_carbon_chain(13, ChainType::Open),
        "tetradec" => Chain::new_carbon_chain(14, ChainType::Open),
        "pentadec" => Chain::new_carbon_chain(15, ChainType::Open),
        "hexadec" => Chain::new_carbon_chain(16, ChainType::Open),
        "heptadec" => Chain::new_carbon_chain(17, ChainType::Open),
        "octadec" => Chain::new_carbon_chain(18, ChainType::Open),
        "nonadec" => Chain::new_carbon_chain(19, ChainType::Open),
        "icos" => Chain::new_carbon_chain(20, ChainType::Open),

        other => {
            todo!("Unimplemented or invalid IUPAC root word {}", other)
        }
    }
}

/// Returns an alcohol functiona group.
/// Can also be used as (oyxgen of) ketone or ether functional group.
pub fn oxygen() -> Chain {
    Chain::new(Atom::Oxygen)
}

/// Returns an aldehyde functional group.
pub fn aldehyde() -> Chain {
    Chain::new(Atom::Carbon).add_atom(Atom::Oxygen, BondType::Double(false))
}

/// Returns a carboxylic acid functional group.
/// Can also be used as ester functional group by adding more atoms.
pub fn carboxylic_acid() -> Chain {
    Chain::new(Atom::Carbon)
        .add_side_chain(&oxygen(), 1, BondType::Double(false))
        .add_atom(Atom::Oxygen, BondType::Single)
}

/// Returns a peroxide/hydroperoxide functional group.
pub fn peroxide() -> Chain {
    Chain::new(Atom::Oxygen).add_atom(Atom::Oxygen, BondType::Single)
}

/// Returns an amine functional group.
/// Can also be used as an imine functional group.
pub fn amine() -> Chain {
    Chain::new(Atom::Nitrogen)
}

/// Returns an amide functional group.
pub fn amide() -> Chain {
    Chain::new(Atom::Carbon)
        .add_side_chain(&oxygen(), 1, BondType::Double(false))
        .add_atom(Atom::Nitrogen, BondType::Single)
}

/// Returns an amidine functional group.
/// Param r2 represents alkyl attatched to nitrogen with double bond.
pub fn amidine(r2: Option<Chain>) -> Chain {
    if let Some(r2) = r2 {
        Chain::new(Atom::Carbon)
            .add_side_chain(
                &Chain::new(Atom::Nitrogen).add_side_chain(&r2, 1, BondType::Single),
                1,
                BondType::Double(false),
            )
            .add_atom(Atom::Nitrogen, BondType::Single)
    } else {
        Chain::new(Atom::Carbon)
            .add_side_chain(&Chain::new(Atom::Nitrogen), 1, BondType::Double(false))
            .add_atom(Atom::Nitrogen, BondType::Single)
    }
}

/// Returns an imide functional group.
pub fn imide() -> Chain {
    Chain::new(Atom::Carbon)
        .add_side_chain(&oxygen(), 1, BondType::Double(false))
        .add_atom(Atom::Nitrogen, BondType::Single)
        .add_atom(Atom::Carbon, BondType::Single)
        .add_side_chain(&oxygen(), 3, BondType::Double(false))
}

/// Returns an azide functional group.
pub fn azide() -> Chain {
    Chain::new(Atom::Nitrogen)
        .add_atom(Atom::Nitrogen, BondType::Double(false))
        .add_atom(Atom::Nitrogen, BondType::Double(false))
}

/// Returns an azo functional group.
pub fn azo() -> Chain {
    Chain::new(Atom::Nitrogen).add_atom(Atom::Nitrogen, BondType::Double(false))
}

/// Returns a cyanate functional group.
pub fn cyanate() -> Chain {
    Chain::new(Atom::Oxygen)
        .add_atom(Atom::Carbon, BondType::Single)
        .add_atom(Atom::Nitrogen, BondType::Triple)
}

/// Returns a isocyanate functional group.
pub fn isocyanate() -> Chain {
    Chain::new(Atom::Nitrogen)
        .add_atom(Atom::Carbon, BondType::Double(false))
        .add_atom(Atom::Oxygen, BondType::Double(false))
}

/// Returns a nitrate functional group.
pub fn nitrate() -> Chain {
    Chain::new(Atom::Oxygen)
        .add_atom(Atom::Nitrogen, BondType::Single)
        .add_side_chain(&oxygen(), 2, BondType::Double(false))
        .add_atom(Atom::Oxygen, BondType::Single)
}

/// Returns a nitrile functional group.
pub fn nitrile() -> Chain {
    Chain::new(Atom::Carbon).add_atom(Atom::Nitrogen, BondType::Triple)
}

/// Returns a isonitrile functional group.
pub fn isonitrile() -> Chain {
    Chain::new(Atom::Nitrogen).add_atom(Atom::Carbon, BondType::Triple)
}

/// Returns a nitrite functional group.
pub fn nitrite() -> Chain {
    Chain::new(Atom::Oxygen)
        .add_atom(Atom::Nitrogen, BondType::Single)
        .add_atom(Atom::Oxygen, BondType::Double(false))
}

/// Returns a nitro functional group.
pub fn nitro() -> Chain {
    Chain::new(Atom::Nitrogen)
        .add_side_chain(&oxygen(), 1, BondType::Double(false))
        .add_atom(Atom::Oxygen, BondType::Single)
}

/// Returns a nitroso functional group.
pub fn nitroso() -> Chain {
    Chain::new(Atom::Nitrogen).add_atom(Atom::Oxygen, BondType::Double(false))
}

/// Returns an oxime functional group.
pub fn oxime() -> Chain {
    Chain::new(Atom::Carbon).add_side_chain(
        &Chain::new(Atom::Nitrogen).add_atom(Atom::Oxygen, BondType::Single),
        1,
        BondType::Double(false),
    )
}

/// Returns a pyridine ring.
pub fn pyridine() -> Chain {
    Chain::new(Atom::Nitrogen)
        .add_atoms(Atom::Carbon, BondType::Double(false), 2)
        .add_atoms(Atom::Carbon, BondType::Double(false), 2)
        .add_atom(Atom::Carbon, BondType::Double(false))
        .convert_chain_type(ChainType::MonoCyclic)
}

/// Returns a carbamate group.
pub fn carbamate() -> Chain {
    Chain::new(Atom::Oxygen)
        .add_atom(Atom::Carbon, BondType::Single)
        .add_side_chain(&oxygen(), 2, BondType::Double(false))
        .add_atom(Atom::Nitrogen, BondType::Single)
}

/// Returns a benzene/phenyl ring.
pub fn benzene() -> Chain {
    Chain::new(Atom::Carbon)
        .add_atom(Atom::Carbon, BondType::Double(false))
        .add_atom(Atom::Carbon, BondType::Single)
        .add_atom(Atom::Carbon, BondType::Double(false))
        .add_atom(Atom::Carbon, BondType::Single)
        .add_atom(Atom::Carbon, BondType::Double(false))
        .convert_chain_type(ChainType::MonoCyclic)
}

/// Returns a toluene ring.
pub fn toluene() -> Chain {
    benzene().add_side_chain(&root_word_chain("meth"), 1, BondType::Single)
}

/// Returns a benzyl ring.
pub fn benzyl() -> Chain {
    root_word_chain("meth").add_side_chain(&benzene(), 1, BondType::Single)
}

/// Returns a phenol ring.
pub fn phenol() -> Chain {
    benzene().add_side_chain(&oxygen(), 1, BondType::Single)
}
