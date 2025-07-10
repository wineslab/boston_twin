from xml.etree import ElementTree as ET
import xml.dom.minidom as minidom
from typing import List, Dict, Any, Union
import numpy as np

def create_mitsuba_xml(
    output_file: str,
    shapes: List[Dict[str, Any]],
    materials: List[Dict[str, Union[str, float]]],
    legacy: bool = False,
) -> str:
    """
    Create a Mitsuba compatible XML file.
    
    Parameters:
    -----------
    output_file : str
        Path to save the XML file
    shapes : list of dict
        List of shape dictionaries with properties like id, filename, material_id
    materials : list of dict
        List of material dictionaries with properties like id, type, thickness
    legacy : bool, optional
        If True, creates XML in legacy format, by default False
        
    Notes:
    ------
    The to_world transformation should be provided as a dictionary with 'rotations' and 'translate' keys:
    {
        'rotations': [{'axis': 'x', 'angle': 120.357}, ...], 
        'translate': [1.287, -5.899, 3.536]
    }
    """
    
    # Create the root element
    root = ET.Element("scene")
    root.set("version", "2.1.0")

    if legacy:
        # Add default comment
        defaults_comment = ET.Comment(
            " Defaults, these can be set via the command line: -Darg=value ")
        root.append(defaults_comment)

        # Add camera comment
        camera_comment = ET.Comment(" Camera and Rendering Parameters ")
        root.append(camera_comment)

        # Add integrator
        integrator = ET.SubElement(root, "integrator")
        integrator.set("type", "path")
        integrator.set("id", "integrator")
        integrator.set("name", "integrator")

        # Add materials comment
        materials_comment = ET.Comment(" Materials ")
        root.append(materials_comment)

        # Add materials in legacy format
        for material in materials:
            bsdf = ET.SubElement(root, "bsdf")
            bsdf.set("type", "twosided")
            mat_id = f"mat-itu_{material['type']}"
            bsdf.set("id", mat_id)
            bsdf.set("name", mat_id)

            nested_bsdf = ET.SubElement(bsdf, "bsdf")
            nested_bsdf.set("type", "diffuse")
            nested_bsdf.set("name", "bsdf")

            rgb = ET.SubElement(nested_bsdf, "rgb")
            rgb.set("value", "0.750000 0.250000 0.500000")
            rgb.set("name", "reflectance")

        # Add emitters comment
        emitters_comment = ET.Comment(" Emitters ")
        root.append(emitters_comment)

        # Add emitter
        emitter = ET.SubElement(root, "emitter")
        emitter.set("type", "constant")
        emitter.set("id", "light")
        emitter.set("name", "light")

        # Add shapes comment
        shapes_comment = ET.Comment(" Shapes ")
        root.append(shapes_comment)

        # Add shapes in legacy format
        for shape in shapes:
            shape_elem = ET.SubElement(root, "shape")
            shape_elem.set("type", "ply")
            shape_elem.set("id", shape["id"])
            shape_elem.set("name", shape["id"])

            # In legacy format, we always add a transform
            transform_elem = ET.SubElement(shape_elem, "transform")
            transform_elem.set("name", "to_world")

            matrix_elem = ET.SubElement(transform_elem, "matrix")
            
            # Handle to_world transformation - convert from dictionary if needed
            if "to_world" in shape:
                if isinstance(shape["to_world"], dict):
                    # Convert dictionary to flattened matrix
                    to_world = transform_dict_to_matrix(shape["to_world"])
                else:
                    # Already a flattened matrix
                    to_world = shape["to_world"]
            else:
                # Default identity matrix
                to_world = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]
                
            # Apply the special legacy matrix ordering
            if len(to_world) == 16:
                # Reorder the matrix elements to match legacy format
                legacy_matrix = [
                    to_world[0], to_world[1], to_world[2], to_world[3],   # First row
                    to_world[8], to_world[10], to_world[9], to_world[11], # Second row (reordered)
                    to_world[4], to_world[6], to_world[5], to_world[7],   # Third row (reordered)
                    to_world[12], to_world[13], to_world[14], to_world[15] # Fourth row
                ]
                matrix_value = " ".join([f"{x:.6f}" for x in legacy_matrix])
            else:
                # Fallback for invalid matrix
                matrix_value = "1.000000 0.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 0.000000 1.000000 0.000000 0.000000 0.000000 0.000000 1.000000"
            
            matrix_elem.set("value", matrix_value)

            filename_elem = ET.SubElement(shape_elem, "string")
            filename_elem.set("name", "filename")
            filename_elem.set("value", shape["filename"])

            normals_elem = ET.SubElement(shape_elem, "boolean")
            normals_elem.set("name", "face_normals")
            normals_elem.set("value", "true")

            ref_elem = ET.SubElement(shape_elem, "ref")
            mat_id = f"mat-itu_{shape['material_id']}"
            ref_elem.set("id", mat_id)
            ref_elem.set("name", "bsdf")

        # Add volumes comment
        volumes_comment = ET.Comment(" Volumes ")
        root.append(volumes_comment)

    else:
        # New format implementation
        materials_comment = ET.Comment(" Materials ")
        root.append(materials_comment)

        # Add materials
        for material in materials:
            bsdf = ET.SubElement(root, "bsdf")
            bsdf.set("type", "itu-radio-material")
            bsdf.set("id", material["id"])

            type_elem = ET.SubElement(bsdf, "string")
            type_elem.set("name", "type")
            type_elem.set("value", material["type"])

            thickness_elem = ET.SubElement(bsdf, "float")
            thickness_elem.set("name", "thickness")
            thickness_elem.set("value", str(material["thickness"]))

        # Add shapes comment
        shapes_comment = ET.Comment(" Shapes ")
        root.append(shapes_comment)

        # Add shapes
        for shape in shapes:
            shape_elem = ET.SubElement(root, "shape")
            shape_elem.set("type", "ply")
            shape_elem.set("id", shape["id"])

            filename_elem = ET.SubElement(shape_elem, "string")
            filename_elem.set("name", "filename")
            filename_elem.set("value", shape["filename"])

            normals_elem = ET.SubElement(shape_elem, "boolean")
            normals_elem.set("name", "face_normals")
            normals_elem.set("value", "true")

            ref_elem = ET.SubElement(shape_elem, "ref")
            ref_elem.set("id", shape["material_id"])
            ref_elem.set("name", "bsdf")

            # Handle to_world transformation
            if "to_world" in shape:
                to_world_elem = ET.SubElement(shape_elem, "transform")
                to_world_elem.set("name", "to_world")

                # Process the transformation dictionary
                if isinstance(shape["to_world"], dict):
                    # Add rotations
                    if "rotations" in shape["to_world"]:
                        for rotation in shape["to_world"]["rotations"]:
                            rotate_elem = ET.SubElement(to_world_elem, "rotate")
                            axis = rotation["axis"]
                            rotate_elem.set(axis, "1")
                            angle = rotation["angle"]
                            rotate_elem.set("angle", f"{angle}")

                    # Add translation
                    if "translate" in shape["to_world"]:
                        translate_elem = ET.SubElement(to_world_elem, "translate")
                        translate_values = shape["to_world"]["translate"]
                        translate_str = " ".join([f"{x:.6f}" for x in translate_values])
                        translate_elem.set("value", translate_str)
                else:
                    # Backward compatibility for flattened matrix
                    matrix_elem = ET.SubElement(to_world_elem, "matrix")
                    matrix_value = " ".join([f"{x:.6f}" for x in shape["to_world"]])
                    matrix_elem.set("value", matrix_value)

    # Convert to a string with pretty formatting
    rough_string = ET.tostring(root, encoding='utf-8')
    reparsed = minidom.parseString(rough_string)

    # Set the indentation style based on format
    indent = "\t" if legacy else "  "
    pretty_xml = reparsed.toprettyxml(indent=indent)

    # Remove extra blank lines that minidom adds
    pretty_lines = [line for line in pretty_xml.split('\n') if line.strip()]
    clean_xml = '\n'.join(pretty_lines)

    # Write to file
    with open(output_file, 'w') as f:
        f.write(clean_xml)

    print(f"Created Mitsuba XML file: {output_file}")
    return clean_xml

def get_transformation_matrix(translation: List[float]=[0,0,0], scale: Union[float,List[float]]=1.0)->List[float]:
    """
    Get the transformation matrix for translation and scaling.
    
    Parameters:
    -----------
    translation : list of float, optional
        List of 3 floats for translation along x, y, z axes, by default [0,0,0]
    scale : float or list of float, optional
        Scaling factor for uniform scaling or list of 3 floats for non-uniform scaling, by default 1.0
    
    Returns:
    --------
    list of float
        List of 16 floats representing the transformation matrix
    """
    if isinstance(scale, float) or isinstance(scale, int):
        scale = [scale, scale, scale]
    return [
        scale[0], 0, 0, translation[0],
        0, scale[1], 0, translation[1],
        0, 0, scale[2], translation[2],
        0, 0, 0, 1
    ]

# Helper function to convert dictionary transformation to 4x4 matrix
def transform_dict_to_matrix(transform_dict):
    # Start with identity matrix
    matrix = np.identity(4)
    
    # Apply rotations if present
    if 'rotations' in transform_dict:
        for rotation in transform_dict['rotations']:
            axis = rotation['axis']
            angle_rad = float(rotation['angle']) * np.pi / 180.0
            
            # Create rotation matrix
            if axis == 'x':
                rot_matrix = np.array([
                    [1, 0, 0, 0],
                    [0, np.cos(angle_rad), -np.sin(angle_rad), 0],
                    [0, np.sin(angle_rad), np.cos(angle_rad), 0],
                    [0, 0, 0, 1]
                ])
            elif axis == 'y':
                rot_matrix = np.array([
                    [np.cos(angle_rad), 0, np.sin(angle_rad), 0],
                    [0, 1, 0, 0],
                    [-np.sin(angle_rad), 0, np.cos(angle_rad), 0],
                    [0, 0, 0, 1]
                ])
            elif axis == 'z':
                rot_matrix = np.array([
                    [np.cos(angle_rad), -np.sin(angle_rad), 0, 0],
                    [np.sin(angle_rad), np.cos(angle_rad), 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]
                ])
            
            # Apply rotation
            matrix = np.matmul(matrix, rot_matrix)
    
    # Apply translation if present
    if 'translate' in transform_dict:
        tx, ty, tz = transform_dict['translate']
        translation_matrix = np.array([
            [1, 0, 0, tx],
            [0, 1, 0, ty],
            [0, 0, 1, tz],
            [0, 0, 0, 1]
        ])
        matrix = np.matmul(matrix, translation_matrix)
        
    # Convert matrix to flattened list
    return matrix.flatten().tolist()