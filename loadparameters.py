from __future__ import absolute_import

import collections
import os
import pprint
import xlrd


class LoadParameters:
    """Class to compute interferograms.
    """

    def __init__(self, sky_spreadsheet='Sky.xlsx', sky_sheet='1point',
      instrument_spreadsheet='FIInS_Instrument_cor3.xlsx'):
        self.result = collections.OrderedDict()
        self.result['substages'] = collections.OrderedDict()

        # access the relevant sheet of the spreadsheet
        book = xlrd.open_workbook(sky_spreadsheet)
        sheet = book.sheet_by_name(sky_sheet)

        sheet_dict = collections.OrderedDict()

        # Sheet has names in col 0, values in columns 1 to n
        # Ignore row 0.

        # read column names
        for row in range(1,sheet.nrows):
            colname = sheet.cell(row, 0)
            sheet_dict[colname.value] = {}
            for col in range(1,sheet.ncols):
                sheet_dict[colname.value][col] = sheet.cell(row,col).value

        self.result['substages']['Sky'] = sheet_dict

        # the instrument spreadsheet
        instrument_book = xlrd.open_workbook(instrument_spreadsheet)
        for sheet in instrument_book.sheets():

            if sheet.name == 'FTSpectrograph':
                sheet_dict = collections.OrderedDict()

                # ignore first row of sheet
                # second row of sheet has column names
                name_row = 1
                # subsequent rows have entries, only one of which is
                # 'selected'
                select_col = None

                # read column names
                for col in range(sheet.ncols):
                    sheet_dict[sheet.cell(name_row,col).value] = \
                      collections.OrderedDict()
                    if 'Selection' in sheet.cell(name_row,col).value:
                        select_col = col

                # sanity check 
                if select_col is None:
                    raise Exception, \
                      'No Selection column in FTSpectrograph spreadsheet'
        
                # read entries from first 'selected' row
                for row in range(name_row+1,sheet.nrows):
                    if int(sheet.cell(row,select_col).value) > 0:
                        for col in range(sheet.ncols):
                            sheet_dict[sheet.cell(name_row,col).value][row] = \
                              sheet.cell(row,col).value
                        break

                # resulting dict should contain the entries from the
                # selected spreadsheet row

                self.result['substages'][sheet.name] = sheet_dict

            elif sheet.name == 'FTSMechanical':
                sheet_dict = collections.OrderedDict()

                # Sheet has names in col 0, values in col 1.
                # Ignore row 0.
                name_col = 0
                entry_col = 1

                # read column names
                for row in range(1,sheet.nrows):
                    sheet_dict[sheet.cell(row,name_col).value] = \
                      collections.OrderedDict()
                    sheet_dict[sheet.cell(row,name_col).value][entry_col] = \
                      sheet.cell(row,entry_col).value

                self.result['substages'][sheet.name] = sheet_dict

            elif sheet.name == 'Telescope':
                sheet_dict = collections.OrderedDict()

                # Sheet has names in col 0, values in col 1, units in col 3.
                # Ignore rows 0 and 1.
                name_col = 0
                entry_col = 1

                # read column names
                for row in range(2,sheet.nrows):
                    sheet_dict[sheet.cell(row,name_col).value] = \
                      sheet.cell(row,entry_col).value

                self.result['substages'][sheet.name] = sheet_dict

            elif sheet.name == 'Interferometer':
                sheet_dict = collections.OrderedDict()

                # ignore first 4 rows of sheet
                # 5th row of sheet has column names
                name_row = 4
                # subsequent rows have entries, only one of which is
                # 'selected'
                select_col = None

                # read column names
                for col in range(sheet.ncols):
                    sheet_dict[sheet.cell(name_row,col).value] = \
                      collections.OrderedDict()
                    if 'Select' in sheet.cell(name_row,col).value:
                        select_col = col

                # sanity check 
                if select_col is None:
                    raise Exception, \
                      'No Selection column in Interferometer spreadsheet'
        
                # read entries from first 'selected' row
                for row in range(name_row+1,sheet.nrows):
                    select_cell = str(sheet.cell(row,select_col).value)
                    select_val = 0
                    try:
                        select_val = int(float(select_cell))
                    except:
                        pass
           
                    if select_val == 1:
                        for col in range(sheet.ncols):
                            sheet_dict[sheet.cell(name_row,col).value][row] = \
                              sheet.cell(row,col).value
                        break

                # resulting dict should contain the entries from the
                # selected spreadsheet row

                self.result['substages'][sheet.name] = sheet_dict

            elif sheet.name == 'ColdOptics':
                sheet_dict = collections.OrderedDict()

                # Sheet has names in col 0, values in col 1.
                # Ignore rows 0 and 1 and final row.
                # Ignore cols 2 and higher - duplicate information.
                name_col = 0
                entry_col = 1

                # read column names
                for row in range(2,sheet.nrows-1):
                    sheet_dict[sheet.cell(row,name_col).value] = \
                      collections.OrderedDict()
                    sheet_dict[sheet.cell(row,name_col).value][entry_col] = \
                      sheet.cell(row,entry_col).value

                self.result['substages'][sheet.name] = sheet_dict

            elif sheet.name == 'Background':
                sheet_dict = collections.OrderedDict()

                # Sheet has names in col 0, values in col 1.
                # Ignore rows 0 and 1.
                name_col = 0
                entry_col = 1

                for row in range(2,sheet.nrows):
                    sheet_dict[sheet.cell(row,name_col).value] = \
                      collections.OrderedDict()
                    sheet_dict[sheet.cell(row,name_col).value][entry_col] = \
                      sheet.cell(row,entry_col).value

                self.result['substages'][sheet.name] = sheet_dict

            elif sheet.name == 'WarmOptics':
                sheet_dict = collections.OrderedDict()

                # Sheet has names in col 0, values in col 1.
                # Ignore rows 0 and 1.
                name_col = 0
                entry_col = 1

                for row in range(2,sheet.nrows):
                    sheet_dict[sheet.cell(row,name_col).value] = \
                      collections.OrderedDict()
                    sheet_dict[sheet.cell(row,name_col).value][entry_col] = \
                      sheet.cell(row,entry_col).value

                self.result['substages'][sheet.name] = sheet_dict

            elif sheet.name == 'Detectors':
                sheet_dict = collections.OrderedDict()

                # Sheet has names in col 0, values in col 1.
                # Ignore rows 0 and 1.
                name_col = 0
                entry_col = 1

                for row in range(2,sheet.nrows):
                    sheet_dict[sheet.cell(row,name_col).value] = \
                      collections.OrderedDict()
                    sheet_dict[sheet.cell(row,name_col).value][entry_col] = \
                      sheet.cell(row,entry_col).value

                self.result['substages'][sheet.name] = sheet_dict

            elif sheet.name == 'SimulatorControl':
                sheet_dict = collections.OrderedDict()

                # Sheet has names in col 0, values in col 1.
                # Ignore row 0.
                name_col = 0
                entry_col = 1

                for row in range(1,sheet.nrows):
                    sheet_dict[sheet.cell(row,name_col).value] = \
                      collections.OrderedDict()
                    sheet_dict[sheet.cell(row,name_col).value][entry_col] = \
                      sheet.cell(row,entry_col).value

                self.result['substages'][sheet.name] = sheet_dict

            else:
                print 'loadparameters skipping', sheet.name

    def run(self):
        return self.result
