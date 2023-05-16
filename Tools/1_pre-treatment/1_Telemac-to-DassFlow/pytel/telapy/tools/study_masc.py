from telapy.api.masc import Mascaret


class StudyMascaret(object):
    def __init__(self, study_files, print_key=0):
        # creation of the mascaret instance
        self.print_key = print_key
        self.masc = Mascaret()
        self.masc.create_mascaret(self.print_key)
        self.study_files = study_files

    def import_model_and_initialisation(self):
        file_type = []
        file_name = []
        for key_val in self.study_files.items():
            if type(key_val[1]) == list:
                for _, sub in enumerate(key_val[1]):
                    file_name.append(sub)
                    file_type.append(key_val[0])
            else:
                file_name.append(key_val[1])
                file_type.append(key_val[0])
        print file_name, len(file_name)
        self.masc.import_model(file_name, file_type)
        self.masc.init_hydro_from_file(self.study_files['lig'])

    def error_message(self):
        self.masc.error_message()

    def run_model(self):

        dt = self.masc.get_double('Model.DT', 0, 0, 0)
        tini = self.masc.get_double('Model.InitTime', 0, 0, 0)
        tfin = self.masc.get_double('Model.MaxCompTime', 0, 0, 0)

        self.masc.compute(tini, tfin, dt)

    def run_model_boucle(self):
        dt = self.masc.get_double('Model.DT', 0, 0, 0)
        tini = self.masc.get_double('Model.InitTime', 0, 0, 0)
        tfin = self.masc.get_double('Model.MaxCompTime', 0, 0, 0)
        t0 = tini
        tplus = t0+dt
        while tplus <= tfin:
            self.masc.compute(t0, tplus, dt)
            t0 = tplus
            tplus = tplus+dt

    def finalize(self):
        del self.masc


def main():
    """
    main function containing an example of a run
    """
    study_files = {'xcas': 'mascaret0.xcas',
                   'geo': 'mascaret0.geo',
                   'res': 'mascaret0.opt',
                   'res_casier': 'mascaret0.opt_casier',
                   'res_liaison': 'mascaret0.opt_liaison',
                   'listing': 'mascaret0.lis',
                   'listing_casier': 'mascaret0.lis_casier',
                   'lig': 'mascaret0.lig',
                   'casier': 'mascaret0.casier',
                   'rep': 'mascaret0_ecr.rep',
                   'listing_liaison': 'mascaret0.lis_liaison',
                   'loi': ['mascaret0_0.loi', 'mascaret0_1.loi']}
    print_key = 0  # 1 pour avoir le listing
    masc = StudyMascaret(study_files, print_key)
    masc.error_message()
    masc.import_model_and_initialisation()
    masc.error_message()
    masc.run_model_boucle()
    masc.error_message()
    masc.finalize()


if __name__ == "__main__":
    main()
