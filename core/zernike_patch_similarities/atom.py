#!/usr/bin/python

import copy
import math

class Atom:
  """Contains attributes for data contained in ATOM tags from PDB-formatted
  files. Additionally, it is able to deal with modifications to fields in
  a non-standard way, like the ones introduced by mark_sur.
  PDB fields: http://www.wwpdb.org/documentation/format33/sect9.html#ATOM
  mark_sur fields: http://northstar-www.dartmouth.edu/doc/insightII/zdockpro/b_fileformats.html
  """
  # Constants for different formattings of the lines
  INPUT_EMPTY = 0
  # INPUT_PDB_BASIC only triggers the load of the basic fields
  INPUT_PDB_BASIC = 1
  INPUT_PDB_STANDARD = 2
  INPUT_MARK_SUR = 3

  # Distance constants (in angstroms)
  INTERFACE_DISTANCE = 5.0
  CONTACT_DISTANCE = 4.0

  # All data members are initialized to None so that
  # we should actually test when using them if they have some sensible value

  # We keep the original line provided to the init method
  # for exporting convenience. Note that that direct calls to the
  # input adding methods will not change this value
  original_line = None

  ##### PDB fields

  # PDB fields 0-indexed. All fields are left-filled with spaces
  serial_number = None # 6-10
  name = None # 12-15
  alternate_location = None # 16
  residue_name = None # 17-19
  # 20 is not specified
  chain_id = None # 21
  residue_sequence_number = None # 22-25
  residue_insertion_code = None # 26
  # 27-29 not specified
  # X,Y,Z coordinates have length 8, 3 decimals strictly
  # and left-filled with spaces
  x = None # 30-37
  y = None # 38-45
  z = None # 46-53
  # The next two are length six with two strict decimals
  occupancy = None # 54-59
  temp_factor = None # 60-65
  # 66-75 not specified
  element_symbol = None # 76-77
  charge = None # 78-79

  ##### mark_sur fields
  
  # The indices are sort of guessed from the output since no
  # reliable documentation has been found.
  mark_sur_atom_type = None # 55-56 (0-unknown, 19-blocked, the rest unicharmm)
  mark_sur_is_exposed = None # 62 (1-0 flag where 1 is exposed)
  mark_sur_vdw_radius = None # 64-67 (length 4, 2 decimals)
  mark_sur_code = None # (?? haven't seen it) PDB code or segment ID string
  mark_sur_charge = None # 77-81 (length five: [sign][digit].[2 decimals]

  # Fields related to additional computations performed on atoms
  is_interface = None
  
  def __init__(self, input_line=None, input_mode=INPUT_EMPTY):
    """The normal use case is for the read process to have the ATOM line
    so it's easier to have a constructor to break it up. Given the variability
    of modifications that are found, from adding extra data at the end to
    violating chain ID rules, the process of reading is broken up into
    sections and the input_mode simply indicates what should be read.
    """
    if (input_mode != Atom.INPUT_EMPTY) and (input_line is not None):
      self.original_line = input_line
      # Assume every format should at least contain these
      self.add_basic_data(input_line)
      if (input_mode == Atom.INPUT_PDB_STANDARD) :
        self.add_non_basic_pdb_data(input_line)
      elif (input_mode == Atom.INPUT_MARK_SUR):
        self.add_mark_sur_data(input_line)

  def add_basic_data(self, input_line):
    """These include fields up to x,y,z coordinates. After that
    point they can have different meanings, like what is done
    by the mark_sur program.
    """
    self.serial_number = self._parse_field(input_line, 6, 5, int)
    self.name = self._parse_field(input_line, 12, 4, str)
    self.alternate_location = self._parse_field(input_line, 16, 1, str)
    self.residue_name = self._parse_field(input_line, 17, 3, str)
    self.chain_id = self._parse_field(input_line, 21, 1, str)
    self.residue_sequence_number = self._parse_field(input_line, 22, 4, int)
    self.residue_insertion_code = self._parse_field(input_line, 26, 1, str)
    self.x = self._parse_field(input_line, 30, 8, float)
    self.y = self._parse_field(input_line, 38, 8, float)
    self.z = self._parse_field(input_line, 46, 8, float)

  def add_non_basic_pdb_data(self, input_line):
    """All PDB fields after x,y,z coordinates
    """
    self.occupancy = self._parse_field(input_line, 54, 6, float)
    self.temp_factor = self._parse_field(input_line, 60, 6, float)
    self.element_symbol = self._parse_field(input_line, 76, 2, str)
    self.charge = self._parse_field(input_line, 78, 2, str)

  def add_mark_sur_data(self, input_line):
    """Parse all instance variables names mark_sur
    """
    self.mark_sur_atom_type = self._parse_field(input_line, 55, 2, str)
    self.mark_sur_is_exposed = self._parse_field(input_line, 62, 1, int)
    self.mark_sur_vdw_radius = self._parse_field(input_line, 64, 4, float)
    # mark_sur_code is not seen in mark_sur output
    self.mark_sur_charge = self._parse_field(input_line, 77, 5, float)

  def _parse_field(self, input_line, start_index, length, type_conversion):
    """Reads the substring and casts to type_conversion. We normally expect
    int, float or str. Returns None if the substring is empty.
    """
    line_section = input_line[start_index : start_index+length].strip()
    if len(line_section) == 0:
      return None
    else:
      return type_conversion(line_section)

  def is_surface(self):
    """True if mark_sur_is_exposed was defined and it's set to 1
    """
    return (self.mark_sur_is_exposed is not None) and self.mark_sur_is_exposed

  def is_closer_than(self, other_atom, threshold):
    """Compute the euclidean distance between self and other atom and,
    if the value is less than threshold, return true.
    """
    # To avoid computing square root, we actually compare the
    # squared distance to the squared threshold
    deltax = self.x - other_atom.x
    deltay = self.y - other_atom.y
    deltaz = self.z - other_atom.z
    return (deltax * deltax + deltay * deltay + deltaz * deltaz) \
           < threshold * threshold

  def shift(self, x, y, z):
    """Shift the atomic coordinates by adding x, y and z to each
    of the corresponding axes.
    """
    self.x += x
    self.y += y
    self.z += z

  def mirror_x(self):
    """This method assumes that the atomic coordinates have been shifted to
    the origin with the sole purpose of generating a mirror image by flipping
    x coordinates to their opposite sign.
    """
    self.x = -self.x

  def set_interface_from_collection(self, atom_collection):
    """An atom is on the interface if its closer than Atom.INTERFACE_DISTANCE
    to any of the atoms in the AtomCollection provided.
    This atom's x,y,z coordinates are compared to every atom in atom_collection
    and the self.is_interface boolean value is set accordingly.
    """
    distance_comparisons = \
        map(lambda other_atom: \
                self.is_closer_than(other_atom, Atom.INTERFACE_DISTANCE),\
            atom_collection.atoms)
    if (reduce(lambda x, y: x or y, distance_comparisons)):
      self.is_interface = True


class AtomCollection:
  """In general this class implements methods that apply to collections
  of atoms, from reading/writing from/to PDB files, to selecting subsets of
  atoms based on their properties.
  """
  atoms = []
  filename = None
  # Only important when each chain is read as a separate atom collection
  chain_id = None

  def __init__(self, input_filename=None, input_mode=Atom.INPUT_PDB_STANDARD):
    """Either create an empty collection if no parameters are provided
    or load atoms from a PDB file (filename). By default we assume it's
    a standard PDB file but the mode can be provided too (see Atom class).
    """
    if (input_filename is not None):
      self.atoms = []
      self.filename = input_filename
      atoms_file = open(input_filename, "r")
      for input_line in atoms_file:
        if input_line.startswith("ATOM"):
          self.atoms.append(Atom(input_line, input_mode))
      atoms_file.close()

  @classmethod
  def new_from_atom_list(cls, atom_list, input_filename=None, chain=None):
    """Factory that creates a new AtomCollection instance from an existing
    list of Atom instances.
    """
    new_atom_collection = cls()
    new_atom_collection.atoms = atom_list
    new_atom_collection.filename = input_filename
    new_atom_collection.chain_id = chain

    return new_atom_collection

  @classmethod
  def new_collections_by_chain(cls, filename, input_mode=Atom.INPUT_PDB_STANDARD):
    """Analog to the main constructor but in this case it will return
    a separate AtomCollection per chain. So, this actually returns a list
    instead of a single instance
    """
    atom_collections = []
    new_atoms = []
    current_chain = ""
    if (input_filename is not None):
      atoms_file = open(input_filename, "r")
      for input_line in atoms_file:
        if input_line.startswith("ATOM"):
          current_atom = Atom(input_line, input_mode)
          if current_chain is "" or current_atom.chain_id == current_chain:
            new_atoms.append(current_atom)
          else:
            atom_collections.append(AtomCollection.new_from_atom_list(\
                new_atoms, filename, current_chain))
            new_atoms = [current_atom]
          current_chain = current_atom.chain_id
      atoms_file.close()
    # insert the last outstanding collection
    atom_collections.append(AtomCollection.new_from_atom_list(\
        new_atoms, filename, current_chain))
    return atom_collections

  @classmethod
  def calculate_interfaces_from_files(cls, filenames, input_mode):
    """Return N AtomCollection instances (N=len(filenames))
    where each one of them contains only the interface atoms.
    Since any number of PDB files are expected, each collection is
    compared against all others and an atom is considered an interface
    one if it's interacting with at least one other PDB.
    """
    atom_collections = []
    for filename in filenames:
      atom_collections.append(cls(filename, input_mode))
    # Mark all interface atoms in all collections
    for base_index in range(len(filenames)):
      for compared_to_index in range(len(filenames)):
        if (base_index != compared_to_index):
          atom_collections[base_index].set_atoms_in_interface(\
              atom_collections[compared_to_index])
    # Interfaces are marked, return just those subsets of atoms
    return map(lambda collection: collection.get_interface_atoms(),\
               atom_collections)

  def length(self):
    """Number of atoms in this collection.
    """
    return len(self.atoms)

  def get_atom_by_serial_number(self, serial_number):
    for atom in self.atoms:
      if atom.serial_number == serial_number:
        return atom
    # If it doesn't find anything
    return None

  def calculate_centroid_distance(self, other_collection):
    this_centroid = self.get_geometric_centroid()
    other_centroid = other_collection.get_geometric_centroid()

    deltax = this_centroid[0] - other_centroid[0]
    deltay = this_centroid[1] - other_centroid[1]
    deltaz = this_centroid[2] - other_centroid[2]
    return math.sqrt(deltax * deltax + deltay * deltay + deltaz * deltaz)

  def get_geometric_centroid(self):
    """Returns a list of 3 numbers that represent the
    x, y, z coordinates of this collection's centroid.
    """
    x = y = z = 0
    for atom in self.atoms:
      x += atom.x
      y += atom.y
      z += atom.z
    return [float(x) / len(self.atoms),\
            float(y) / len(self.atoms),\
            float(z) / len(self.atoms)] 

  def get_atoms_shifted_to_origin(self):
    min_x = min_y = min_z = float("inf")
    for atom in self.atoms:
      if atom.x < min_x:
        min_x =  atom.x
      if atom.y < min_y:
        min_y =  atom.y
      if atom.z < min_z:
        min_z =  atom.z
    shifted_atoms = copy.deepcopy(self)
    map(lambda atom : atom.shift(-min_x, -min_y, -min_z), shifted_atoms.atoms)
    return shifted_atoms
#      print str(atom.x) + " " + str(atom.y) + " " + str(atom.z)
#      print str(shifted_atom.x) + " " + str(shifted_atom.y) + " " + str(shifted_atom.z)
#    return AtomCollection.new_from_atom_list(shifted_atoms);

  def get_x_mirrored_atoms(self):
    mirrored_atoms = copy.deepcopy(self)
    map(lambda atom: atom.mirror_x(), mirrored_atoms.atoms)
    return mirrored_atoms

  def get_surface_atoms(self):
    """A new AtomCollection is returned with all the atoms whose
    is_surface member is set to true.
    """
    return AtomCollection.new_from_atom_list(\
        filter(lambda atom: atom.is_surface(), self.atoms))

  def get_interface_atoms(self):
    """A new AtomCollection is returned with all the atoms whose
    is_interface member is set to true.
    """
    return AtomCollection.new_from_atom_list(\
        filter(lambda atom : atom.is_interface, self.atoms))

  def get_atoms_close_to(self, other_atom, threshold):
    """Compare all atoms against other atoms. Return a list
    with the ones that are closer than threshold.
    """
    close_atoms = AtomCollection.new_from_atom_list(\
        filter(lambda atom: atom.is_closer_than(other_atom, threshold),\
        self.atoms))
    # Update the filename based on the other atom used.
    # This serves for exporting purposes.
    filename_parts = self.filename.split(".", 1)
    close_atoms.filename = "%s-atom%05dresidue%04d.%s" % \
        (filename_parts[0],other_atom.serial_number,\
        other_atom.residue_sequence_number,filename_parts[1])
    return close_atoms

  def get_surface_patches(self, threshold):
    """An atom collection annotated by mark_sur contains surface
    atoms s_1, s_2,...., s_n. For each s_i we create a new AtomCollection
    that contains s_i itself and all other atoms closer than threshold to
    s_i's x,y,z coordinates.
    The return value is thus a list where each position i is the AtomCollection
    that represents the patch for s_i.
    """
    return map(lambda surface_atom:\
                   self.get_atoms_close_to(surface_atom, threshold),\
                   self.get_surface_atoms().atoms)

  def set_atoms_in_interface(self, atom_collection):
    """Set the is_interface attribute for all self.atoms by comparing
    each of them to the atom_collection provided.
    """
    for atom in self.atoms:
      atom.set_interface_from_collection(atom_collection)

  def count_overlap(self, other_collection):
    """Returns the number of atoms that are by both this instance
    and other_collection.
    """
    self_serials = set(map(lambda atom: atom.serial_number, self.atoms))
    other_serials = set(map(lambda atom: atom.serial_number,\
                            other_collection.atoms))
    # this calculates the intersection
    common_serials = self_serials & other_serials
    return len(common_serials)

  def write_to_file(self, output_filename=None):
    """Writes each of the atoms to the output_filename using
    the original_line attribute from each atom. If output_filename
    is None, the filename in this AtomCollection is used.
    The only values that are taken from the instance variables are the
    x, y, z coordinates. This is due to the variability of usage that
    some fields have in different PDB modified formats (like mark_sur,
    alternate position conventions, etc).
    """
    if output_filename is None:
      output_filename = self.filename
    output_file = open(output_filename, "w");
    for atom in self.atoms:
      output_line = "%s%8.3f%8.3f%8.3f%s" % \
          (atom.original_line[:30], atom.x, atom.y, atom.z,\
           atom.original_line[54:])
      output_file.write(output_line)
    output_file.close()

  def print_all(self):
    """Print all original lines in the atom collection.
    """
    print self.filename
    for atom in self.atoms:
      print atom.original_line,
