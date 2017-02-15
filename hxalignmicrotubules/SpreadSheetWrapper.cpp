#include <hxalignmicrotubules/SpreadSheetWrapper.h>

HX_INIT_CLASS(SpreadSheetWrapper, HxSpreadSheet);

const char* SpreadSheetWrapper::ROW_NAME = "rowName\0";

SpreadSheetWrapper::SpreadSheetWrapper() : mRowHashMap(-1) {}

SpreadSheetWrapper::~SpreadSheetWrapper(void) {}

void SpreadSheetWrapper::addEntry(const char* tableID, const char* columnID,
                                  const char* rowID, const char* value) {
    McHandle<HxSpreadSheet> theSpreadSheet = getSpreadSheet();

    const int tableIdx = getTableIdx(tableID);

    int rowIdx = getRowIdx(tableID, columnID, rowID);

    const int columnIdx = getColumnIdx(tableID, columnID);

    HxSpreadSheet::Column* actColumn =
        theSpreadSheet->column(columnIdx, tableIdx);

    actColumn->setValue(rowIdx, value);
}

void SpreadSheetWrapper::addEntry(const char* tableID, const char* columnID,
                                  const char* rowID, float value) {
    McHandle<HxSpreadSheet> theSpreadSheet = getSpreadSheet();

    const int tableIdx = getTableIdx(tableID);

    int rowIdx = getRowIdx(tableID, columnID, rowID);

    const int columnIdx =
        getColumnIdx(tableID, columnID, HxSpreadSheet::Column::STRING);

    McString valueAsString;
    valueAsString.printf("%f", value);

    HxSpreadSheet::Column* actColumn =
        theSpreadSheet->column(columnIdx, tableIdx);

    actColumn->setValue(rowIdx, valueAsString);
}

void SpreadSheetWrapper::getEntry(const char* tableID, const char* columnID,
                                  const char* rowID, McString& entry) {

    entry.clear();
    McHandle<HxSpreadSheet> theSpreadSheet = getSpreadSheet();

    const int tableIdx = getTableIdx(tableID);

    const int columnIdx = getColumnIdx(tableID, columnID);

    int rowIdx = getRowIdx(tableID, columnID, rowID);

    HxSpreadSheet::Column* actColumn =
        theSpreadSheet->column(columnIdx, tableIdx);

    if (actColumn->type == HxSpreadSheet::Column::FLOAT) {
        entry.printf("%f", actColumn->floatValue(rowIdx));
    } else
        entry = McString(actColumn->stringValue(rowIdx));
}

void SpreadSheetWrapper::initRowColumn(const char* tableID) {
    const int rowColumnIdx =
        getColumnIdx(tableID, McString(ROW_NAME).getString());
    mcassert(rowColumnIdx == 0);
}

int SpreadSheetWrapper::getTableIdx(const char* tableID) {
    McHandle<HxSpreadSheet> theSpreadSheet = getSpreadSheet();
    int tableIdx = theSpreadSheet->findTable(tableID);
    if (tableIdx == -1) {
        theSpreadSheet->addTable(tableID);
        tableIdx = theSpreadSheet->findTable(tableID);
        initRowColumn(tableID);
    }
    return tableIdx;
}

int SpreadSheetWrapper::parse(Tcl_Interp* t, int argc, char** argv) {
    char* cmd = argv[1];
    if (CMD("setAValue")) {
        ASSERTARG(6);
        addEntry(argv[2], argv[3], argv[4], argv[5]);
        return TCL_OK;
    }
    if (CMD("getAValue")) {
        ASSERTARG(5);
        McString entry;
        getEntry(argv[2], argv[3], argv[4], entry);
        Tcl_VaSetResult(t, "%s", entry.getString());
        return TCL_OK;
    } else {
        return HxSpreadSheet::parse(t, argc, argv);
    }
}

int SpreadSheetWrapper::getColumnIdx(const char* tableID, const char* columnID,
                                     int type) {

    McHandle<HxSpreadSheet> theSpreadSheet = getSpreadSheet();
    int tableIdx = getTableIdx(tableID);
    int columnIdx = theSpreadSheet->findColumn(columnID, type, tableIdx);
    if (columnIdx == -1) {
        theSpreadSheet->addColumn(columnID, type, tableIdx);
        columnIdx = theSpreadSheet->findColumn(columnID, type, tableIdx);
    }
    return columnIdx;
}

int SpreadSheetWrapper::getRowIdx(const char* tableID, const char* columnID,
                                  const char* rowID) {

    McHandle<HxSpreadSheet> theSpreadSheet = getSpreadSheet();
    int tableIdx = getTableIdx(tableID);

    McString rowKey = createRowKey(tableID, rowID);

    RowIdx* rowIdx = mRowHashMap.lookup(rowKey);
    if (rowIdx == 0) {
        int numRows = theSpreadSheet->nRows(tableIdx);
        theSpreadSheet->setNumRows(numRows + 1, tableIdx);
        mRowHashMap.insert(rowKey, RowIdx(numRows));
        addEntry(tableID, McString(ROW_NAME).getString(), rowID, rowID);
        rowIdx = mRowHashMap.lookup(rowKey);
    }

    mcassert(rowIdx);
    return rowIdx->mIndex;
}

McString SpreadSheetWrapper::createRowKey(const char* tableID,
                                          const char* rowID) {
    McString rowKey(rowID);
    rowKey.remove(rowKey.size() - 1, 1);
    McString mcTableID(tableID);
    rowKey.append(mcTableID.size(), tableID);
    return rowKey;
}

McHandle<HxSpreadSheet> SpreadSheetWrapper::getSpreadSheet() { return this; }

int SpreadSheetWrapper::getNumRows(const char* tableID) {
    McHandle<HxSpreadSheet> theSpreadSheet = getSpreadSheet();
    const int tableIdx = getTableIdx(tableID);
    if (tableIdx >= 0)
        return theSpreadSheet->nRows(tableIdx);
    else
        return -1;
}
