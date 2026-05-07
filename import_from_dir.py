import sys
import importlib.util
import os

def import_from_directory(directory_path, module_name):
    # Add directory to the system path
    sys.path.append(directory_path)
    
    # Load the module
    module_path = os.path.join(directory_path, f"{module_name}.py")
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    
    return module

# Example usage
# external_module = import_from_directory("/path/to/directory", "external_script")
# external_module.some_function()


