#pragma once

#include <mclib/McHandle.h>
#include <mclib/McHashTable.h>

#ifdef HX_AMIRA5_COMPAT
#include <hxstatistics/HxSpreadSheet.h>
#else
#include <hxspreadsheet/HxSpreadSheet.h>
#endif

#include <hxalignmicrotubules/api.h>

class HXALIGNMICROTUBULES_API SpreadSheetWrapper : public HxSpreadSheet {

    HX_HEADER(SpreadSheetWrapper);

  public:
    class RowIdx {
      public:
        RowIdx(int theIndex = -1) { mIndex = theIndex; };

        int mIndex;
    };

    SpreadSheetWrapper();

    ~SpreadSheetWrapper(void);

    static const char* ROW_NAME;

    void addEntry(const char* tableID, const char* columnId, const char* rowID,
                  const char* value);
    void addEntry(const char* tableID, const char* columnId, const char* rowID,
                  float value);

    void getEntry(const char* tableID, const char* columnID, const char* rowID,
                  McString& value);
    virtual int parse(Tcl_Interp* t, int argc, char** argv);

    int getNumRows(const char* tableID);

  private:
    int getTableIdx(const char* tableID);

    int getColumnIdx(const char* tableID, const char* columnID,
                     int type = HxSpreadSheet::Column::STRING);

    int getRowIdx(const char* tableID, const char* columnID, const char* rowID);

    McHandle<HxSpreadSheet> getSpreadSheet();

    McHashTable<McString, SpreadSheetWrapper::RowIdx> mRowHashMap;

    McString createRowKey(const char* tableID, const char* rowID);

    void initRowColumn(const char* tableID);
};
