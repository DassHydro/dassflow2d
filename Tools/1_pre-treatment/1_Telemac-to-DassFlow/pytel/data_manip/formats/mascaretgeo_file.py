from collections import OrderedDict
import os.path

from data_manip.formats.mascaret_file import Reach, Section
from utils.exceptions import MascaretException


class MascaretGeoFile:
    """
    Parse Mascaret geometry file (geo/geoC/georef/georefC)
    Handles multiple reaches
    TODO: handle major/minor bed

    Attributs:
    - file_name (str) file name
    - fformat (str) file format ('opt' or 'rub')

    - has_ref (bool): has X and Y coordinates for points
        (and hydraulic axis position)
    - has_layers (bool): has sediment layers
    - nlayers (int): number of layers
    - layer_names (list): list of layer names
    """
    OUTPUT_FLOAT_FMT = '%.6f'

    def __init__(self, file_name, fformat=None, mode='read'):
        """
        @param file_name (str) file name
        @param fformat (str) file format ('opt' or 'rub')
        @param mode (str) define the mode for the class,
            'read' by default to read a file,
            anything else to create a file
        """
        self.file_name = file_name
        self.reaches = OrderedDict()

        # Layers for sediments (Courlis)
        self.nlayers = 0
        self.layer_names = []

        if mode == 'read':
            # File format information
            if fformat is None:
                self.fformat = os.path.splitext(file_name)[1][1:]
            else:
                self.fformat = fformat.lower().strip()
            if self.fformat not in ('geo', 'georef'):
                raise NotImplementedError(
                    'Format `%s` not supported,\
                     only geo and georef formats are supported as input' %
                    self.fformat)
            self.has_ref = 'ref' in self.fformat

            # Layers for sediments (Courlis)
            self.has_layers = self.fformat.endswith('C')

            # Load file content
            self.load()

    def load(self):
        """
        Load Mascaret geometry file (geo/georef)
        """
        with open(self.file_name, 'r') as filein:
            reach = None
            reach_id_curr = 1
            section_id = 0
            section_name = ''
            section_pk = -1.0
            dist, x_list, y_list, z_list = [], [], [], []
            xa, ya = None, None

            for line in filein:
                if line.startswith('#'):
                    # Ignore comment line
                    pass
                elif line.upper().startswith('PROFIL'):
                    if dist:
                        # Add previous Section
                        section = Section(section_id, section_pk, section_name)
                        if xa is not None and ya is not None:
                            section.set_axis(xa, ya)
                        if self.has_ref:
                            section.set_points_from_xyz(x_list, y_list, z_list)
                        else:
                            section.set_points_from_trans(dist, z_list)
                        reach.add_section(section)

                    if self.has_ref:
                        _, reach_name, section_name, pk_str, x1, y1, x2, y2,\
                         _, xa, ya = line.split()
                        xa = float(xa)
                        ya = float(ya)
                    else:
                        _, reach_name, section_name, pk_str = line.split()

                    # Create first reach for initialisation
                    if reach is None:
                        reach = Reach(reach_id_curr, reach_name)
                        self.add_reach(reach)
                        reach_id_curr += 1

                    # Create a new reach if necessary
                    if reach_name != reach.name:
                        reach = Reach(reach_id_curr, reach_name)
                        self.add_reach(reach)
                        reach_id_curr += 1

                    # Reset variables to store section
                    section_pk = float(pk_str)
                    dist, x_list, y_list, z_list = [], [], [], []
                    section_id += 1
                else:
                    if self.has_ref:
                        dist_str, z_str, _, x1, y1 = line.split()
                        x_list.append(float(x1))
                        y_list.append(float(y1))
                    else:
                        dist_str, z_str, _ = line.split()

                    # Add new point to current section
                    dist.append(float(dist_str))
                    z_list.append(float(z_str))

            # Add last section
            section = Section(section_id, section_pk, section_name)
            if xa is not None and ya is not None:
                section.set_axis(xa, ya)
            if self.has_ref:
                section.set_points_from_xyz(x_list, y_list, z_list)
            else:
                section.set_points_from_trans(dist, z_list)
            reach.add_section(section)

    def save(self, output_file_name):
        """
        Save Mascaret geometry file (geo/georef)
        @param output_file_name (str) output file name
        """
        fformat = os.path.splitext(output_file_name)[1]
        if fformat == '.geo':
            ref, layers = False, False
        elif fformat == '.georef':
            ref, layers = True, False
        elif fformat == '.geoC':
            ref, layers = False, True
        else:  # georefC
            ref, layers = True, True

        if ref and not self.has_ref:
            raise MascaretException('Could not write `%s` format without\
             any geo-referenced data' % fformat)

        with open(output_file_name, 'w') as fileout:
            for _, reach in self.reaches.items():

                for sec in reach:
                    positions_str = ''
                    if ref:
                        # Get river_banks and `AXE` coordinates if necessary
                        xa, ya = sec.axis
                        positions_str += ' %f %f %f %f' %\
                            (sec.x[0], sec.y[0], sec.x[-1], sec.y[-1])
                        positions_str += ' AXE %f %f' % (xa, ya)

                    # Write profile header
                    fileout.write(
                        'Profil %s %s %f%s\n' %
                        (reach.name, sec.name, sec.pk, positions_str))

                    # Write points and layers if necessary
                    if not ref and not layers:
                        for dist, z in zip(sec.distances, sec.z):
                            fileout.write('%f %f B\n' % (dist, z))

                    elif ref and not layers:
                        for dist, x, y, z in zip(sec.distances,
                                                 sec.x, sec.y, sec.z):
                            fileout.write('%f %f B %f %f\n' % (dist, z, x, y))

                    elif not ref and layers:
                        for i, (dist, z) in enumerate(zip(sec.distances,
                                                          sec.z)):
                            if self.nlayers == 0:
                                layers_str = ''
                            else:
                                layers_str = ' ' +\
                                    ' '.join(
                                        [MascaretGeoFile.OUTPUT_FLOAT_FMT % zl
                                            for zl in sec.layers_elev[:, i]])
                            fileout.write('%f %f%s B\n' %
                                          (dist, z, layers_str))

                    elif ref and layers:
                        for i, (dist, x, y, z) in enumerate(zip(sec.distances,
                                                            sec.x, sec.y,
                                                            sec.z)):
                            if self.nlayers == 0:
                                layers_str = ''
                            else:
                                layers_str = ' ' + ' '\
                                    .join([MascaretGeoFile.OUTPUT_FLOAT_FMT %
                                          zl for zl in sec.layers_elev[:, i]])
                            fileout.write('%f %f%s B %f %f\n' %
                                          (dist, z, layers_str, x, y))

    def __repr__(self):
        return 'MascaretGeoFile: %s' % self.file_name

    def add_reach(self, reach):
        """
        Add a single reach
        @param reach (Reach) reach to add
        """
        self.reaches[reach.id] = reach

    def add_constant_layer(self, name, thickness):
        self.has_layers = True
        self.nlayers += 1
        self.layer_names.append(name)
        for _, reach in self.reaches.items():
            for section in reach:
                thickness_table = [thickness for i in range(section.nb_points)]
                section.add_layer(thickness_table)

    def summary(self):
        txt = '~> %s\n' % self
        for _, reach in self.reaches.items():
            txt += '    - %s\n' % reach
            for section in reach:
                txt += '        %i) %s\n' % (section.id, section)
        return txt


if __name__ == '__main__':
    # Parse every MascaretGeoFile in examples/mascaret and display a summary
    import os
    from utils.files import recursive_glob
    try:
        geo_files = recursive_glob(os.path.join(os.environ['HOMETEL'],
                                   'examples', 'mascaret'), '*.geo')
        geo_files += recursive_glob(os.path.join(os.environ['HOMETEL'],
                                    'examples', 'mascaret'), 'geometrie')
        for file_name in sorted(geo_files):
            geo_file = MascaretGeoFile(file_name, 'geo')
            print(geo_file.summary())
    except MascaretException as e:
        print(str(e))
