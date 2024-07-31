"""
A GUI for the DART project.
"""
import sys
import json
from PyQt5.QtWidgets import QApplication, QMainWindow, QTabWidget, QWidget, QVBoxLayout, QLabel, QLineEdit, QPushButton, QFileDialog, QScrollArea, QComboBox, QFrame, QHBoxLayout, QSpacerItem, QSizePolicy
from PyQt5.QtGui import QIntValidator
from PyQt5.QtCore import Qt

metals = ('Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn')  # as option for the metal center in the GUI

from DARTassembler.src.assembly.Assembly_Input import allowed_topologies, _verbose, _optimization_movie, _concatenate_xyz, _overwrite_output_path, _output_path, _complex_name_length, _name, _input_path, _max_num_complexes, _ligand_choice, _isomers, _optimisation, _random_seed, _total_charge, _topology, _element, _oxidation_state, _complex_name_appendix, _geometry_modifier_filepath, _bidentate_rotator, _same_ligand_keyword

assembler_options = {
    _output_path: {
                'object': 'path',
                'info': 'Path to the output folder.',
                'default': None
    }}
advanced_options = {
    _overwrite_output_path: {
                'object': (True, False),
                'info': 'Whether to overwrite the output if it already exists.',
                'default': False
    },
    _optimization_movie: {
                'object': (True, False),
                'info': 'Whether to output a movie (concatenated .xyz) of the forcefield optimization. Set to false to save disk space.',
                'default': True
    },
    _concatenate_xyz: {
                'object': (True, False),
                'info': 'Whether to save a file with all final structures concatenated.',
                'default': True
    },
    _verbose: {
                'object': (0, 1, 2, 3),
                'info': 'How much output to print. 0: no output, 1: minimal output, 2: normal output, 3: verbose output.',
                'default': 2
    },
    _complex_name_length: {
                'object': (50, 20, 15, 12, 8, 5),
                'info': 'Length of the complex name.',
                'default': 8
    }
}
batch_options = {
    _name: {
                'object': str,
                'info': 'Name of the batch.',
                'default': None
    },
    _topology: {
                'object': tuple(allowed_topologies),
                'info': 'Geometry of the assembled complexes. Currently supported: Octahedral (3-2-1, 4-1-1, 5-1) or square-planar (2-1-1, 2-2). The numbers denote the denticities of each ligand.',
                'default': '3-2-1'
    },
    _input_path: {
                'object': f'Either single path or list of [path, "{_same_ligand_keyword}"].',
                'info': f'Path to the ligand database. Either single path or list of [path, "{_same_ligand_keyword}"].',
                'default': None
    },
    _ligand_choice: {
                'object': ('random', 'all'),
                'info': 'How to choose the ligands.',
                'default': 'random'
    },
    _max_num_complexes: {
                'object': int,
                'info': 'Maximum number of complexes to generate.',
                'default': 100
    },
    _element: {
                'object': metals,
                'info': 'Chemical symbol of the desired metal center.',
                'default': 'Pd'
    },
    _oxidation_state: {
                    'object': (7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3),
                    'info': 'Oxidation state of the desired metal center.',
                    'default': 3,
    },
    _total_charge: {
                    'object': (5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5),
                    'info': 'Total charge of the complex.',
                    'default': 0,
    },
    _optimisation: {
                    'object': (True, False),
                    'info': 'Whether to optimize the structures after generation with a force field.',
                    'default': False,
    },
    _isomers: {
                    'object': ('lowest_energy', 'all'),
                    'info': 'Which isomers to generate.',
                    'default': 'lowest_energy',
    },
    _bidentate_rotator: {
                    'object': ('auto', 'horseshoe', 'slab'),
                    'info': 'How to rotate the bidentate ligands.',
                    'default': 'auto',
    },
    _geometry_modifier_filepath: {
                    'object': 'path',
                    'info': 'Path to the geometry modifier file. If not given, no geometry modification is performed.',
                    'default': None,
    },
    _random_seed: {
                    'object': int,
                    'info': 'Random seed for reproducibility.',
                    'default': 0,
    },
    _complex_name_appendix: {
                    'object': str,
                    'info': 'String to append to the randomly generated complex name.',
                    'default': None
    }
}

ligand_options = {
    'input_ligand_db_path': {
        'object': 'path or string',
        'info': 'Path to the input ligand database or "MetaLig".',
        'default': 'MetaLig'
    },
    'output_ligand_db_path': {
        'object': 'path',
        'info': 'Path to the output ligand database.',
        'default': 'filtered_ligand_db.json'
    }
}
instructions = ('must_contain_and_only_contain', 'must_at_least_contain', 'must_exclude', 'must_only_contain_in_any_amount')
apply_to_denticities = {
            'object': 'list[int] or None',
            'info': 'List of denticities to apply this filter to. If empty, applies to all denticities.',
            'default': None,
        }
ligand_filters = {
    'denticities': {
        'denticities': {
            'object': 'list[int]',
            'info': 'Only keep ligands with these denticities.',
            'default': 'Placeholder - 1, 2, 3',
        },
    },
    'ligand_charges': {
        'ligand_charges': {
            'object': 'list[int]',
            'info': 'Only keep ligands with these charges.',
            'default': 'Placeholder - -1, 0, 1',
        },
        'apply_to_denticities': apply_to_denticities
    },
    'ligand_composition': {
        'elements': {
            'object': 'list[str]',
            'info': 'Elements to apply this filter to.',
            'default': 'Placeholder - C, N, O',
        },
        'instruction': {
            'object': instructions,
            'info': 'Instruction for how to apply this filter.',
            'default': 'must_contain_and_only_contain',
        },
        'apply_to_denticities': apply_to_denticities
    },
    'coordinating_atoms_composition': {
        'elements': {
            'object': 'list[str]',
            'info': 'Elements to apply this filter to.',
            'default': 'Placeholder - C, N, O',
        },
        'instruction': {
            'object': instructions,
            'info': 'Instruction for how to apply this filter.',
            'default': 'must_contain_and_only_contain',
        },
        'apply_to_denticities': apply_to_denticities
    },
    'number_of_atoms': {
        'min': {
            'object': 'int or None',
            'info': 'Minimum number of atoms. If empty, defaults to 0.',
            'default': 'Placeholder - 0',
        },
        'max': {
            'object': 'int or None',
            'info': 'Maximum number of atoms. If empty, defaults to infinity.',
            'default': 'Placeholder - 40',
        },
        'apply_to_denticities': apply_to_denticities
    },
    'molecular_weight': {
        'min': {
            'object': 'float or None',
            'info': 'Minimum molecular weight (in g/mol). If empty, defaults to 0.',
            'default': 'Placeholder - 0',
        },
        'max': {
            'object': 'float or None',
            'info': 'Maximum molecular weight (in g/mol). If empty, defaults to infinity.',
            'default': 'Placeholder - 500',
        },
        'apply_to_denticities': apply_to_denticities
    },
    'metal_donor_bond_lengths': {
        'min': {
            'object': 'float or None',
            'info': 'Minimum bond length (in Angstrom). If empty, defaults to 0.',
            'default': 'Placeholder - 1.2',
        },
        'max': {
            'object': 'float or None',
            'info': 'Maximum bond length (in Angstrom). If empty, defaults to infinity.',
            'default': None,
        },
        'apply_to_denticities': apply_to_denticities
    },
    'interatomic_distances': {
        'min': {
            'object': 'float or None',
            'info': 'Minimum interatomic distance (in Angstrom). If empty, defaults to 0.',
            'default': 'Placeholder - 0.6',
        },
        'max': {
            'object': 'float or None',
            'info': 'Maximum interatomic distance (in Angstrom). If empty, defaults to infinity.',
            'default': None,
        },
        'apply_to_denticities': apply_to_denticities
    },
    'planarity': {
        'min': {
            'object': 'float or None',
            'info': 'Minimum planarity score. If empty, defaults to 0.0.',
            'default': 'Placeholder - 0.0',
        },
        'max': {
            'object': 'float or None',
            'info': 'Maximum planarity score. If empty, defaults to 1.0.',
            'default': 'Placeholder - 1.0',
        },
        'apply_to_denticities': apply_to_denticities
    },
    'occurrences': {
        'min': {
            'object': 'int or None',
            'info': 'Minimum number of occurrences. If empty, defaults to 0.',
            'default': 'Placeholder - 0',
        },
        'max': {
            'object': 'int or None',
            'info': 'Maximum number of occurrences. If empty, defaults to infinity.',
            'default': None,
        },
        'apply_to_denticities': apply_to_denticities
    },
    'metal_ligand_binding_history': {
        'metal_ligand_binding_history': {
            'object': 'list[str]',
            'info': 'List of metals to keep.',
            'default': 'Placeholder - Pd, Pt, Ru',
        },
        'apply_to_denticities': apply_to_denticities
    },
    'remove_ligands_with_adjacent_coordinating_atoms': {
        'remove_ligands_with_adjacent_coordinating_atoms': {
            'object': 'bool',
            'info': 'true or false. If false, filter will be ignored.',
            'default': True,
        }
    },
    'remove_ligands_with_beta_hydrogens': {
        'remove_ligands_with_beta_hydrogens': {
            'object': 'bool',
            'info': 'true or false. If false, filter will be ignored.',
            'default': True,
        }
    },
    'graph_IDs': {
        'graph_IDs': {
            'object': 'list[str]',
            'info': 'List of graph IDs to keep.',
            'default': 'Placeholder - a2b7bbb6ca4ce36dc3147760335e7374, 53b7a3d91a1be6e167a3975bb7921206',
        }
    },
    'remove_ligands_with_missing_bond_orders': {
        'remove_ligands_with_missing_bond_orders': {
            'object': 'bool',
            'info': 'true or false. If false, filter will be ignored.',
            'default': True,
        }
    },
    'atomic_neighbors': {
        'atom': {
            'object': 'str',
            'info': 'Chemical element of the central atom.',
            'default': 'Placeholder - C',
        },
        'neighbors': {
            'object': 'str',
            'info': 'List of chemical elements or stoichiometry.',
            'default': 'Placeholder - H2',
        },
        'apply_to_denticities': apply_to_denticities
    },
    'smarts': {
        'smarts': {
            'object': 'str',
            'info': 'SMARTS pattern to match. Must be enclosed in quotes. Example: "[C&H2]".',
            'default': 'Placeholder - "[C&H2]"',
        },
        'should_contain': {
            'object': (True, False),
            'info': 'If True, the ligand must contain the SMARTS pattern to pass. If False, the ligand must not contain the SMARTS pattern to pass.',
            'default': True,
        },
        'include_metal': {
            'object': (True, False),
            'info': 'If True, the ligand structure will contain the metal center "Cu" connected to the coordinating atoms.',
            'default': True,
        },
        'apply_to_denticities': apply_to_denticities
    }
}



class GUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("DART Assembler - Graphical User Interface")
        self.setGeometry(100, 100, 600, 700)

        self.tabs = QTabWidget()
        self.setCentralWidget(self.tabs)

        self.modules = ["ligandfilters", "assembler", "dbinfo", "concat"]
        self.module_methods = {
            "ligandfilters": self.setup_ligandfilters,
            "assembler": self.setup_assembler,
            "dbinfo": self.setup_dbinfo,
            "concat": self.setup_concat
        }

        for module in self.modules:
            tab = QWidget()
            self.tabs.addTab(tab, module.capitalize())
            self.module_methods[module](tab)

    def browse_file(self, line_edit):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(self, "Open File", "", "All Files (*)", options=options)
        if file_path:
            line_edit.setText(file_path)

    def save_input(self, module_name, data):
        with open(f"{module_name}_input.json", "w") as json_file:
            json.dump(data, json_file, indent=4)
        print(f"Saved {module_name}_input.json")

    def setup_ligandfilters(self, tab):
        layout = QVBoxLayout()

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_content = QWidget()
        scroll_layout = QVBoxLayout(scroll_content)

        self.ligandfilter_fields = {}

        def add_widget(layout, label_text, option):
            container = QHBoxLayout()
            label = QLabel(f"{label_text}:")
            container.addWidget(label)

            info_icon = QLabel("?")
            info_icon.setToolTip(option['info'])
            info_icon.setStyleSheet("color: yellow; font-weight: bold;")
            info_icon.setCursor(Qt.PointingHandCursor)
            container.addWidget(info_icon)

            entry = QLineEdit()
            if option['object'] == 'path':
                browse_button = QPushButton("Browse")
                browse_button.clicked.connect(lambda: self.browse_file(entry))
                container.addWidget(entry)
                container.addWidget(browse_button)
            elif isinstance(option['object'], tuple):
                combobox = QComboBox()
                combobox.addItems([str(item) for item in option['object']])
                container.addWidget(combobox)
                self.ligandfilter_fields[label_text] = combobox
                if option['default'] is not None and not option['default'].startswith('Placeholder - '):
                    combobox.setCurrentText(str(option['default']))
            elif option['object'] == int:
                entry.setValidator(QIntValidator())
                container.addWidget(entry)
            elif 'list' in str(option['object']):
                entry.setPlaceholderText("Comma-separated values")
                container.addWidget(entry)
            else:
                container.addWidget(entry)

            if option['default'] is not None:
                if isinstance(option['object'], tuple):
                    if not option['default'].startswith('Placeholder - '):
                        combobox.setCurrentText(str(option['default']))
                else:
                    if option['default'].startswith('Placeholder - '):
                        entry.setPlaceholderText(option['default'].replace('Placeholder - ', ''))
                    else:
                        entry.setText(str(option['default']))

            self.ligandfilter_fields[label_text] = entry
            layout.addLayout(container)

        # Add general ligandfilter options
        for key, option in ligand_options.items():
            add_widget(scroll_layout, key, option)

        # Add visual separator
        separator1 = QFrame()
        separator1.setFrameShape(QFrame.HLine)
        separator1.setFrameShadow(QFrame.Sunken)
        scroll_layout.addWidget(separator1)

        # Filter section
        filter_label = QLabel("Filters:")
        filter_label.setStyleSheet("font-weight: bold;")
        filter_label.setAlignment(Qt.AlignCenter)
        scroll_layout.addWidget(filter_label)

        self.filter_frames = []

        def update_filter_titles():
            for idx, (filter_frame, filter_type) in enumerate(self.filter_frames, start=1):
                title_label = filter_frame.findChild(QLabel, "filterTitle")
                if title_label:
                    title_label.setText(f"Filter {idx}: {filter_type}")

        def add_filter():
            filter_type = filter_dropdown.currentText()
            filter_options = ligand_filters[filter_type]

            filter_frame = QFrame()
            filter_layout = QVBoxLayout(filter_frame)

            title_container = QHBoxLayout()
            title_label = QLabel(f"Filter {len(self.filter_frames) + 1}: {filter_type}")
            title_label.setObjectName("filterTitle")
            title_label.setStyleSheet("font-weight: bold;")
            title_container.addWidget(title_label)

            remove_button = QPushButton("Remove Filter")
            remove_button.setStyleSheet("font-size: 10px; color: red;")
            remove_button.clicked.connect(lambda: remove_filter(filter_frame))
            title_container.addWidget(remove_button)
            title_container.addStretch()

            filter_layout.addLayout(title_container)

            for key, option in filter_options.items():
                add_widget(filter_layout, key, option)

            separator = QFrame()
            separator.setFrameShape(QFrame.HLine)
            separator.setFrameShadow(QFrame.Sunken)
            filter_layout.addWidget(separator)

            self.filter_frames.append((filter_frame, filter_type))
            scroll_layout.insertWidget(scroll_layout.count() - 2, filter_frame)
            scroll_content.adjustSize()
            update_filter_titles()

        def remove_filter(filter_frame):
            for frame, filter_type in self.filter_frames:
                if frame == filter_frame:
                    self.filter_frames.remove((frame, filter_type))
                    break
            filter_frame.setParent(None)
            filter_frame.deleteLater()
            scroll_content.adjustSize()
            update_filter_titles()

        # Add button to add filters
        bottom_container = QFrame()
        bottom_layout = QHBoxLayout(bottom_container)
        bottom_layout.setContentsMargins(0, 0, 0, 0)

        filter_dropdown = QComboBox()
        filter_dropdown.addItems(ligand_filters.keys())
        bottom_layout.addWidget(filter_dropdown)

        add_filter_button = QPushButton("Add Filter")
        add_filter_button.setStyleSheet("font-weight: bold; font-size: 14px; color: green;")
        add_filter_button.clicked.connect(add_filter)
        bottom_layout.addWidget(add_filter_button)

        save_button = QPushButton("Save")
        save_button.setStyleSheet("font-weight: bold; font-size: 14px; color: blue;")
        save_button.clicked.connect(self.collect_ligandfilters_data)
        bottom_layout.addWidget(save_button)

        scroll_layout.addWidget(bottom_container)

        scroll_area.setWidget(scroll_content)
        layout.addWidget(scroll_area)
        layout.addWidget(bottom_container)
        tab.setLayout(layout)

    def collect_ligandfilters_data(self):
        ligandfilters_data = {field_name: entry.text() for field_name, entry in self.ligandfilter_fields.items() if
                              field_name not in ligand_filters}
        filters = []

        for frame, filter_type in self.filter_frames:
            filter_data = {"type": filter_type}
            filter_options = ligand_filters[filter_type]
            for field_name, entry in self.ligandfilter_fields.items():
                if field_name in filter_options:
                    option_type = filter_options[field_name]['object']
                    if 'list' in str(option_type):
                        if option_type == 'list[int]':
                            filter_data[field_name] = [int(x.strip()) for x in entry.text().split(',') if x.strip()]
                        elif option_type == 'list[str]':
                            filter_data[field_name] = [x.strip() for x in entry.text().split(',') if x.strip()]
                    else:
                        filter_data[field_name] = entry.text()
            filters.append(filter_data)

        ligandfilters_data["filters"] = filters
        self.save_input("ligandfilters", ligandfilters_data)

    def setup_assembler(self, tab):
        layout = QVBoxLayout()

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_content = QWidget()
        scroll_layout = QVBoxLayout(scroll_content)

        self.assembler_fields = {}

        def add_widget(layout, label_text, option):
            container = QHBoxLayout()
            label = QLabel(f"{label_text}:")
            container.addWidget(label)

            info_icon = QLabel("?")
            info_icon.setToolTip(option['info'])
            info_icon.setStyleSheet("color: yellow; font-weight: bold;")
            info_icon.setCursor(Qt.PointingHandCursor)
            container.addWidget(info_icon)

            if option['object'] == 'path':
                entry = QLineEdit()
                browse_button = QPushButton("Browse")
                browse_button.clicked.connect(lambda: self.browse_file(entry))
                container.addWidget(entry)
                container.addWidget(browse_button)
                self.assembler_fields[label_text] = entry
            elif isinstance(option['object'], tuple):
                combobox = QComboBox()
                combobox.addItems([str(item) for item in option['object']])
                container.addWidget(combobox)
                self.assembler_fields[label_text] = combobox
            elif option['object'] == int:
                entry = QLineEdit()
                entry.setValidator(QIntValidator())
                container.addWidget(entry)
                self.assembler_fields[label_text] = entry
            else:
                entry = QLineEdit()
                container.addWidget(entry)
                self.assembler_fields[label_text] = entry

            if option['default'] is not None:
                if isinstance(option['object'], tuple):
                    combobox.setCurrentText(str(option['default']))
                else:
                    entry.setText(str(option['default']))

            layout.addLayout(container)

        def add_ligand_db_path_widget(layout, ligand_paths):
            container = QWidget()
            container_layout = QHBoxLayout(container)
            entry = QLineEdit()
            browse_button = QPushButton("Browse")
            browse_button.clicked.connect(lambda: self.browse_file(entry))
            container_layout.addWidget(entry)
            container_layout.addWidget(browse_button)

            same_ligand_button = QPushButton("Same Ligand")
            same_ligand_button.clicked.connect(lambda: entry.setText('same_ligand_as_previous'))
            container_layout.addWidget(same_ligand_button)

            remove_button = QPushButton("Remove")
            remove_button.clicked.connect(lambda: remove_ligand_db_path_widget(container, layout, ligand_paths))
            container_layout.addWidget(remove_button)

            ligand_paths.append(entry)
            layout.addWidget(container)

        def remove_ligand_db_path_widget(container, ligand_paths_layout, ligand_paths):
            for i, entry in enumerate(ligand_paths):
                if entry.parentWidget() == container:
                    ligand_paths.pop(i)
                    break

            # Remove all child widgets from the container
            while container.layout().count():
                item = container.layout().takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.setParent(None)
                    widget.deleteLater()

            # Remove the container itself
            ligand_paths_layout.removeWidget(container)
            container.setParent(None)
            container.deleteLater()

            # Force layout update
            ligand_paths_layout.update()

        # Add assembler options
        for key, option in assembler_options.items():
            add_widget(scroll_layout, key, option)

        # Add visual separator
        separator1 = QFrame()
        separator1.setFrameShape(QFrame.HLine)
        separator1.setFrameShadow(QFrame.Sunken)
        scroll_layout.addWidget(separator1)

        # Batches section
        batches_label = QLabel("Batches:")
        batches_label.setStyleSheet("font-weight: bold;")
        scroll_layout.addWidget(batches_label)
        self.batch_frames = []

        def update_batch_titles():
            for idx, (batch_frame, _) in enumerate(self.batch_frames, start=1):
                title_label = batch_frame.findChild(QLabel, "batchTitle")
                if title_label:
                    title_label.setText(f"Batch {idx}")

        def add_batch():
            batch_frame = QFrame()
            batch_layout = QVBoxLayout(batch_frame)

            title_label = QLabel(f"Batch {len(self.batch_frames) + 1}")
            title_label.setObjectName("batchTitle")
            title_label.setStyleSheet("font-weight: bold;")
            batch_layout.addWidget(title_label)

            for key, option in batch_options.items():
                if key == 'ligand_db_paths':
                    ligand_paths_layout = QVBoxLayout()
                    add_ligand_db_path_button = QPushButton("Add Ligand DB Path")
                    ligand_paths = []
                    add_ligand_db_path_button.clicked.connect(
                        lambda: add_ligand_db_path_widget(ligand_paths_layout, ligand_paths))
                    ligand_paths_layout.addWidget(add_ligand_db_path_button)
                    batch_layout.addLayout(ligand_paths_layout)
                    self.assembler_fields[key] = ligand_paths
                else:
                    add_widget(batch_layout, key, option)

            separator = QFrame()
            separator.setFrameShape(QFrame.HLine)
            separator.setFrameShadow(QFrame.Sunken)
            batch_layout.addWidget(separator)

            self.batch_frames.append((batch_frame, {key: self.assembler_fields[key] for key in batch_options.keys()}))
            # Insert before the batch buttons container
            scroll_layout.insertWidget(scroll_layout.indexOf(batch_buttons_container), batch_frame)
            scroll_content.adjustSize()
            update_batch_titles()

        def remove_batch():
            if self.batch_frames:
                batch_frame, _ = self.batch_frames.pop()
                batch_frame.setParent(None)
                scroll_content.adjustSize()
                update_batch_titles()

        # Add batch and remove batch buttons
        batch_buttons_layout = QHBoxLayout()
        add_batch_button = QPushButton("Add Batch")
        add_batch_button.clicked.connect(add_batch)
        batch_buttons_layout.addWidget(add_batch_button)

        remove_batch_button = QPushButton("Remove Batch")
        remove_batch_button.clicked.connect(remove_batch)
        batch_buttons_layout.addWidget(remove_batch_button)

        batch_buttons_container = QFrame()
        batch_buttons_container.setLayout(batch_buttons_layout)
        scroll_layout.addWidget(batch_buttons_container)

        # Add visual separator
        separator2 = QFrame()
        separator2.setFrameShape(QFrame.HLine)
        separator2.setFrameShadow(QFrame.Sunken)
        scroll_layout.addWidget(separator2)

        # Advanced options
        advanced_options_container = QFrame()
        advanced_options_layout = QVBoxLayout(advanced_options_container)
        advanced_options_container.setVisible(False)

        for key, option in advanced_options.items():
            add_widget(advanced_options_layout, key, option)

        advanced_options_toggle_button = QPushButton("Show Advanced Options")

        def toggle_advanced_options():
            if advanced_options_container.isVisible():
                advanced_options_container.setVisible(False)
                advanced_options_toggle_button.setText("Show Advanced Options")
            else:
                advanced_options_container.setVisible(True)
                advanced_options_toggle_button.setText("Hide Advanced Options")

        advanced_options_toggle_button.clicked.connect(toggle_advanced_options)
        scroll_layout.addWidget(advanced_options_toggle_button)
        scroll_layout.addWidget(advanced_options_container)

        save_button = QPushButton("Save")
        save_button.clicked.connect(self.collect_assembler_data)
        scroll_layout.addWidget(save_button)

        scroll_area.setWidget(scroll_content)
        layout.addWidget(scroll_area)
        tab.setLayout(layout)

    def collect_assembler_data(self):
        assembler_data = {field_name: entry.text() for field_name, entry in self.assembler_fields.items() if
                          field_name != 'ligand_db_paths'}
        batches = []

        for frame, batch_field_entries in self.batch_frames:
            batch_data = {}
            for field_name, entry in batch_field_entries.items():
                if field_name == 'ligand_db_paths':
                    batch_data[field_name] = [e.text() for e in entry]
                else:
                    batch_data[field_name] = entry.text()
            batches.append(batch_data)

        assembler_data["batches"] = batches
        self.save_input("assembler", assembler_data)

    def setup_dbinfo(self, tab):
        layout = QVBoxLayout()

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_content = QWidget()
        scroll_layout = QVBoxLayout(scroll_content)

        label = QLabel("Dbinfo Input Path:")
        scroll_layout.addWidget(label)

        input_entry = QLineEdit()
        scroll_layout.addWidget(input_entry)

        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(lambda: self.browse_file(input_entry))
        scroll_layout.addWidget(browse_button)

        save_button = QPushButton("Save")
        save_button.clicked.connect(lambda: self.save_input("dbinfo", {"input_path": input_entry.text()}))
        scroll_layout.addWidget(save_button)

        scroll_area.setWidget(scroll_content)
        layout.addWidget(scroll_area)
        tab.setLayout(layout)

    def setup_concat(self, tab):
        layout = QVBoxLayout()

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_content = QWidget()
        scroll_layout = QVBoxLayout(scroll_content)

        label = QLabel("Concat Input Path:")
        scroll_layout.addWidget(label)

        input_entry = QLineEdit()
        scroll_layout.addWidget(input_entry)

        browse_button = QPushButton("Browse")
        browse_button.clicked.connect(lambda: self.browse_file(input_entry))
        scroll_layout.addWidget(browse_button)

        save_button = QPushButton("Save")
        save_button.clicked.connect(lambda: self.save_input("concat", {"input_path": input_entry.text()}))
        scroll_layout.addWidget(save_button)

        scroll_area.setWidget(scroll_content)
        layout.addWidget(scroll_area)
        tab.setLayout(layout)


# Main function to start the GUI
def main():
    app = QApplication(sys.argv)
    gui = GUI()
    gui.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
