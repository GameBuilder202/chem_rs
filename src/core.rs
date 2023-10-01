//! This module has the barebones functionality for
//! building molecules and their chemfig and SMILES.

/// Any struct implementing this trait can export their
/// data to chemfig code.
///
/// Note: This should only return the chemfig code without
/// the \chemfig{...} around it.
pub trait ToChemfig {
    /// Self explanatory, returns the chemfig representation.
    fn chemfig_repr(&self) -> String;
}

/// Any struct implementing this trait can export their
/// data to SMILES strings.
pub trait ToSMILES {
    /// Self explanatory, returns the SMILES string.
    fn smiles(&self) -> String;
}

/// This enum basically represents the periodic table of this library.
#[allow(missing_docs)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Atom {
    Carbon,
    Nitrogen,
    Oxygen,
    Sulphur,
}

impl Atom {
    fn get_valency(&self, bonds: u64) -> u64 {
        match self {
            Self::Carbon => 4,
            Self::Nitrogen => 3,
            Self::Oxygen => 2,
            Self::Sulphur => {
                if bonds <= 2 {
                    2
                } else if bonds <= 4 {
                    4
                } else {
                    6
                }
            }
        }
    }

    fn get_attached_hydrogens(&self, bonds: u64) -> u64 {
        let valency = self.get_valency(bonds);
        if bonds > valency {
            panic!("Too many bonds to {:?} atom", self)
        }
        valency - bonds
    }
}

impl ToChemfig for Atom {
    fn chemfig_repr(&self) -> String {
        match self {
            Self::Carbon => String::new(),
            Self::Nitrogen => String::from("N"),
            Self::Oxygen => String::from("O"),
            Self::Sulphur => String::from("S"),
        }
    }
}

impl ToSMILES for Atom {
    fn smiles(&self) -> String {
        match self {
            Self::Carbon => String::from("C"),
            Self::Nitrogen => String::from("N"),
            Self::Oxygen => String::from("O"),
            Self::Sulphur => String::from("S"),
        }
    }
}

/// This enum represents the type of the chain and how it is laid out.
///
/// The arguments for [`ChainType::PolyCyclic`] and [`ChainType::Spiro`] directly
/// correspond to the IUPAC nomenclature for these compounds
#[derive(Debug, Clone)]
pub enum ChainType {
    /// Represents an open chain.
    Open,

    /// Represents closed chains with 1 ring.
    MonoCyclic,

    /// Represents polycyclic chains which have more than 1 ring, such as bicyclo, tricyclo etc.
    ///
    /// The arguments directly correspond to the IUPAC nomenclature.
    ///
    /// Note: See macro [`polycyclic`] if you want to create this variant with an easier syntax.
    ///
    /// # Examples
    /// ```
    /// use chem_rs::core::ChainType;
    ///
    /// // bicyclo[3.2.0]heptane
    /// let chain_type1 = ChainType::PolyCyclic(3, 2, 0, vec![]);
    ///
    /// // tetracyclo[6.6.3.0(2,7).0(9,14)]heptadecane
    /// let chain_type2 = ChainType::PolyCyclic(6, 6, 3, vec![(0, 2, 7), (0, 9, 14)]);
    /// ```
    PolyCyclic(u64, u64, u64, Vec<(u64, u64, u64)>),

    /// Represents a spiro compound. The arguments directly
    /// correspond to the IUPAC nomenclature.
    ///
    /// # Examples
    /// ```
    /// use chem_rs::core::ChainType;
    ///
    /// // spiro[4.3]octane
    /// let chain_type1 = ChainType::Spiro(4, 3);
    ///
    /// // spiro[5.5]undecane
    /// let chain_type2 = ChainType::Spiro(5, 5);
    /// ```
    Spiro(u64, u64),
}

/// A tiny macro to make creation of [`ChainType::PolyCyclic`] easier
#[macro_export]
macro_rules! polycyclic {
    // [$major:literal . $minor:literal . $bridge_len:literal] => {{
    // 	use chem_rs::core::ChainType;
    // 	ChainType::PolyCyclic($major, $minor, $bridge_len, vec![])
    // }};
    [$major:literal . $minor:literal . $bridge_len:literal $(. $len:literal($start:literal, $end:literal))*] => {{
        use chem_rs::core::ChainType;
        ChainType::PolyCyclic($major, $minor, $bridge_len, vec![
            $(($len, $start, $end)),+
        ])
    }};
}

/// Represents how many bonds an atom of the chain has
/// to the *next* atom of the chain or to a side chain.
///
/// *Note: Side chains in chem_rs do not always correspond
/// to alkyl substituenta.*
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum BondType {
    /// Reserved for the last atom of an open chain or of
    /// non-cyclic side chains.
    None,

    /// Represents a single bond to the next atom in the chain
    /// or to a side chain.
    Single,

    /// Represents a double bond to the next atom in the chain
    /// or to a side chain. The boolean argument corresponds to
    /// whether this double bond is a cis/trans configuration,
    /// where `false` represents a cis double bond.
    Double(bool),

    /// Represents a triple bond to the next atom in the chain.
    Triple,
}

impl BondType {
    fn num_bonds(&self) -> u64 {
        match self {
            Self::None => 0,
            Self::Single => 1,
            Self::Double(_) => 2,
            Self::Triple => 3,
        }
    }
}

impl ToChemfig for BondType {
    fn chemfig_repr(&self) -> String {
        match self {
            Self::None => String::from(""),
            Self::Single => String::from("-"),
            Self::Double(_) => String::from("="),
            Self::Triple => String::from("~"),
        }
    }
}

impl ToSMILES for BondType {
    fn smiles(&self) -> String {
        match self {
            Self::None => String::from(""),
            Self::Single => String::from(""),
            Self::Double(_) => String::from("="),
            Self::Triple => String::from("#"),
        }
    }
}

/// Represents a chain of atoms, however it can also be
/// just 1 atom. Usually this will be a chain of carbon
/// atoms, however there is no restriction on what atoms
/// you can chain.
#[derive(Clone)]
pub struct Chain {
    atoms: Vec<(Atom, BondType)>,
    side_chains: Vec<Vec<(Chain, BondType)>>,

    chain_type: ChainType,
}

impl Chain {
    // Creating new chains

    /// Creates a new chain with the given atom as the starting atom.
    ///
    /// # Examples
    /// ```
    /// use chem_rs::core::{Chain, Atom};
    ///
    /// let methane = Chain::new(Atom::Carbon);
    /// let ammonia = Chain::new(Atom::Nitrogen);
    /// ```
    pub fn new(start_atom: Atom) -> Self {
        Self {
            atoms: vec![(start_atom, BondType::None)],
            side_chains: vec![vec![]],
            chain_type: ChainType::Open,
        }
    }

    /// Creates a new carbon chain with the specified length and chain type.
    ///
    /// The polycyclic or spiro configuration's length must match the input length.
    ///
    /// # Examples
    /// ```
    /// use chem_rs::core::{Chain, ChainType};
    ///
    /// let propane = Chain::new_carbon_chain(3, ChainType::Open);
    /// let cyclohexane = Chain::new_carbon_chain(6, ChainType::MonoCyclic);
    /// ```
    pub fn new_carbon_chain(len: u64, chain_type: ChainType) -> Self {
        let chain = match chain_type.clone() {
            ChainType::Open => [
                vec![(Atom::Carbon, BondType::Single); (len - 1) as usize],
                vec![(Atom::Carbon, BondType::None)],
            ]
            .concat(),

            ChainType::MonoCyclic => vec![(Atom::Carbon, BondType::Single); len as usize],

            ChainType::PolyCyclic(major, minor, bridge_len, bridges) => {
                Self::check_polycyclic_config(len, major, minor, bridge_len, &bridges);
                vec![(Atom::Carbon, BondType::Single); len as usize]
            }

            ChainType::Spiro(minor, major) => {
                Self::check_spiro_config(len, minor, major);
                vec![(Atom::Carbon, BondType::Single); len as usize]
            }
        };

        Self {
            atoms: chain,
            side_chains: vec![vec![]; len as usize],
            chain_type,
        }
    }

    // Functions for modifying the chain

    /// Adds the given atom to the chain with the input bond type
    /// being the number of bonds to the previous atom. This function
    /// can be chained with other chain modifying functions to build a
    /// chain.
    ///
    /// # Examples
    /// ```
    /// use chem_rs::core::{Chain, Atom, BondType};
    ///
    /// let hydrogen_cyanide = Chain::new(Atom::Carbon)
    ///                         .add_atom(Atom::Nitrogen, BondType::Triple);
    /// let propene = Chain::new(Atom::Carbon)
    ///                .add_atom(Atom::Carbon, BondType::Double(false))
    ///                .add_atom(Atom::Carbon, BondType::Single);
    /// ```
    pub fn add_atom(mut self, atom: Atom, bond_type: BondType) -> Self {
        let len = self.atoms.len();
        self.atoms[len - 1].1 = bond_type;
        self.atoms.push((atom, BondType::None));
        self.side_chains.push(vec![]);

        self
    }

    /// Adds the given atom n times to the chain with the input bond type
    /// being the number of bonds to the first atom added. This function can
    /// be chained with other chain modifying functions to build a chain.
    ///
    /// # Examples
    /// ```
    /// use chem_rs::core::{Chain, Atom, ChainType, BondType};
    ///
    /// let dodeca_4_en_7_yne = Chain::new(Atom::Carbon)
    ///     .add_atoms(Atom::Carbon, BondType::Single, 3)
    ///     .add_atoms(Atom::Carbon, BondType::Double, 3)
    ///     .add_atoms(Atom::Carbon, BondType::Triple, 4);
    /// ```
    pub fn add_atoms(self, atom: Atom, bond_type: BondType, n: usize) -> Self {
        let mut res = self.add_atom(atom, bond_type);
        for _ in 0..(n - 1) {
            res = res.add_atom(atom, BondType::Single)
        }
        res
    }

    /// Adds the given chain as a side chain to this chain at the specified
    /// locant with the number of bonds to the chain being the `bond_type` argument.
    /// This function can be chained with other chain modifying functions to build a chain.
    ///
    /// # Examples
    /// ```
    /// use chem_rs::core::{Chain, Atom, ChainType, BondType};
    ///
    /// let isobutane = Chain::new_carbon_chain(3, ChainType::Open)
    ///                  .add_side_chain(
    ///                     &Chain::new(Atom::Carbon), 2, BondType::Single
    ///                  );
    /// // This shows an example of a side chain not being an alkyl group.
    /// let acetone = Chain::new_carbon_chain(3, ChainType::Open)
    ///                .add_side_chain(
    ///                     &Chain::new(Atom::Oxygen), 2, BondType::Double(false)
    ///                );
    /// ```
    pub fn add_side_chain(mut self, chain: &Chain, locant: u64, bond_type: BondType) -> Self {
        self.side_chains[(locant - 1) as usize].push((chain.clone(), bond_type));

        self
    }

    /// Converts the chain type of this chain, doing necessary checks beforehand
    pub fn convert_chain_type(mut self, chain_type: ChainType) -> Self {
        match &chain_type {
            ChainType::Open => {}
            ChainType::MonoCyclic => {
                let len = self.atoms.len();
                self.atoms[len - 1].1 = BondType::Single
            }

            ChainType::PolyCyclic(major, minor, bridge_len, bridges) => {
                Self::check_polycyclic_config(
                    self.atoms.len() as u64,
                    *major,
                    *minor,
                    *bridge_len,
                    bridges,
                )
            }

            ChainType::Spiro(minor, major) => {
                Self::check_spiro_config(self.atoms.len() as u64, *minor, *major)
            }
        }

        self.chain_type = chain_type;
        self
    }

    fn check_polycyclic_config(
        len: u64,
        major: u64,
        minor: u64,
        bridge_len: u64,
        bridges: &Vec<(u64, u64, u64)>,
    ) {
        let cycle_len = major
            + minor
            + bridge_len
            + bridges.iter().map(|bridge_info| bridge_info.0).sum::<u64>()
            + 2;
        assert_eq!(
            len,
            cycle_len,
            "Polycyclic configuration does not match input length!\nPolycyclic configuration: {}\nInput length: {}",
            format!("[{}.{}.{}", major, minor, bridge_len)
                + if !bridges.is_empty() {
                    ".".to_owned() + &bridges.iter().map(|bridge| format!("{}({},{})", bridge.0, bridge.1, bridge.2))
                        .collect::<Vec<String>>()
                        .join(".")
                } else { String::new() }.as_str() + "]",
            len
        );
    }

    fn check_spiro_config(len: u64, minor: u64, major: u64) {
        let cycle_len = minor + major + 1;
        assert_eq!(
            len,
            cycle_len,
            "Spiro configuration does not match input length!\nSpiro configuration: {}\nInput length: {}",
            format_args!("[{}.{}]", minor, major),
            len
        );
    }

    fn draw_chemfig(&self, angle_offset: i64, prev_bond: BondType) -> String {
        let mut s = String::new();
        match self.chain_type {
            ChainType::Open => {
                let mut angle = angle_offset + 30;
                let mut sign = -1i64;
                for i in 0..self.atoms.len() {
                    let atom = &self.atoms[i];
                    let side_chains = &self.side_chains[i];

                    s += &atom.0.chemfig_repr();

                    let hydrogens = atom.0.get_attached_hydrogens(
                        side_chains.iter().map(|(_, b)| b.num_bonds()).sum::<u64>()
                            + if i > 0 {
                                self.atoms[i - 1].1.num_bonds()
                            } else {
                                0
                            }
                            + self.atoms[i].1.num_bonds()
                            + prev_bond.num_bonds(),
                    );
                    if hydrogens > 0 && atom.0 != Atom::Carbon {
                        if hydrogens > 1 {
                            s += format!("H_{}", hydrogens).as_str()
                        } else {
                            s += "H"
                        }
                    }

                    if atom.1 == BondType::None {
                        break;
                    }

                    for side_chain in side_chains {
                        s += format!(
                            "([:{},,1]{}{})",
                            angle + 120 * sign,
                            side_chain.1.chemfig_repr(),
                            side_chain
                                .0
                                .draw_chemfig(angle + 120 * sign - 90, side_chain.1)
                        )
                        .as_str()
                    }

                    s += format!("{}[:{},,1]", atom.1.chemfig_repr(), angle).as_str();

                    angle += 60 * sign;
                    sign *= -1;
                }
            }

            ChainType::MonoCyclic => {
                let len = self.atoms.len();
                let inc = (360 / len) as i64;
                let mut angle = -(90 + inc);
                s += format!(
                    "{}{}*{{{}}}(-",
                    {
                        let hydrogens = self.atoms[0].0.get_attached_hydrogens(
                            self.side_chains[0]
                                .iter()
                                .map(|(_, b)| b.num_bonds())
                                .sum::<u64>()
                                + self.atoms[0].1.num_bonds()
                                + self.atoms[self.atoms.len() - 1].1.num_bonds()
                                + prev_bond.num_bonds(),
                        );

                        if hydrogens > 0 && self.atoms[0].0 != Atom::Carbon {
                            if hydrogens > 1 {
                                format!("H_{}", hydrogens)
                            } else {
                                String::from("H")
                            }
                        } else {
                            String::new()
                        }
                    },
                    self.atoms[0].0.chemfig_repr(),
                    self.atoms.len()
                )
                .as_str();
                for i in 1..self.atoms.len() {
                    let atom = &self.atoms[i];
                    let side_chains = &self.side_chains[i];
                    s += atom.0.chemfig_repr().as_str();

                    let hydrogens = atom.0.get_attached_hydrogens(
                        side_chains.iter().map(|(_, b)| b.num_bonds()).sum::<u64>()
                            + if i > 0 {
                                self.atoms[i - 1].1.num_bonds()
                            } else {
                                0
                            }
                            + self.atoms[i].1.num_bonds()
                            + prev_bond.num_bonds(),
                    );
                    if hydrogens > 0 && atom.0 != Atom::Carbon {
                        if hydrogens > 1 {
                            s += format!("H_{}", hydrogens).as_str()
                        } else {
                            s += "H"
                        }
                    }

                    for side_chain in side_chains {
                        s += format!(
                            "({}[,,1]{})",
                            side_chain.1.chemfig_repr(),
                            side_chain
                                .0
                                .draw_chemfig(angle - 90 + angle_offset, side_chain.1)
                        )
                        .as_str();
                    }

                    s += format!("{}[,,1]", atom.1.chemfig_repr()).as_str();

                    angle += inc
                }
                s += ")";
            }

            _ => todo!(),
        }

        s
    }

    fn get_smiles(&self, prev_bond: BondType) -> String {
        let mut s = String::new();
        match self.chain_type {
            ChainType::Open => {
                let mut prev_bond = prev_bond;
                for i in 0..self.atoms.len() {
                    let atom = &self.atoms[i];
                    let side_chains = &self.side_chains[i];

                    let valency = atom
                        .0
                        .get_attached_hydrogens(prev_bond.num_bonds() + atom.1.num_bonds());
                    if valency > 0 && atom.0 != Atom::Carbon {
                        s += format!("[{}H{}]", atom.0.smiles(), valency).as_str();
                    } else {
                        s += &atom.0.smiles();
                    }
                    for side_chain in side_chains {
                        s += format!(
                            "({}{})",
                            side_chain.1.smiles(),
                            side_chain.0.get_smiles(side_chain.1)
                        )
                        .as_str()
                    }
                    s += atom.1.smiles().as_str();

                    prev_bond = atom.1
                }
            }

            ChainType::MonoCyclic => {
                let first = &self.atoms[0];
                s += (first.0.smiles() + "1").as_str();
                for i in 1..self.atoms.len() {
                    let atom = &self.atoms[i];
                    let side_chains = &self.side_chains[i];

                    s += &atom.0.smiles();
                    for side_chain in side_chains {
                        s +=
                            format!("({}{})", side_chain.1.smiles(), side_chain.0.smiles()).as_str()
                    }
                    s += atom.1.smiles().as_str()
                }
                s += "1"
            }

            _ => todo!(),
        }

        s
    }
}

impl ToChemfig for Chain {
    /// Returns the chemfig representation of the molecule.
    ///
    /// # Examples
    /// ```
    /// use chem_rs::core::{Chain, Atom, ChainType, BondType};
    ///
    /// let cyclohexane = Chain::new_carbon_chain(6, ChainType::MonoCyclic);
    /// println!("cyclohexane chemfig: {}", cyclohexane.chemfig_repr());
    ///
    /// let hydrogen_cyanide = Chain::new(Atom::Carbon).add_atom(Atom::Nitrogen, BondType::Triple);
    /// println!("hydrogen_cyanide chemfig: {}", hydrogen_cyanide.chemfig_repr());
    /// ```
    #[inline]
    fn chemfig_repr(&self) -> String {
        self.draw_chemfig(0, BondType::None)
    }
}

impl ToSMILES for Chain {
    /// Returns the SMILES of the molecule.
    ///
    /// # Examples
    /// ```
    /// use chem_rs::core::{Chain, Atom, ChainType, BondType};
    ///
    /// let cyclohexane = Chain::new_carbon_chain(6, ChainType::MonoCyclic);
    /// println!("cyclohexane SMILES: {}", cyclohexane.smiles());
    ///
    /// let hydrogen_cyanide = Chain::new(Atom::Carbon).add_atom(Atom::Nitrogen, BondType::Triple);
    /// println!("hydrogen_cyanide SMILES: {}", hydrogen_cyanide.smiles());
    /// ```
    fn smiles(&self) -> String {
        self.get_smiles(BondType::None)
    }
}
