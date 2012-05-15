
import os

from nesoni import grace, io, config, reference_directory

class Working(io.Workspace):
    def set_reference(self, path):
        self.update_param(reference=self.path_as_relative_path(path))

    def setup_reference(self, filenames):
        if len(filenames) == 1 and os.path.isdir(filenames[0]):
            self.set_reference(filenames[0])
            return
        
        path = self.object_filename('reference')
        reference_directory.Make_reference(path, filenames).run()
        self.set_reference(path)

    def get_reference(self):
        if 'reference' in self.param:
            path = self.relative_path_as_path(self.param['reference'])
        else:
            path = self.working_dir
        
        return reference_directory.Reference(path, must_exist=True)



