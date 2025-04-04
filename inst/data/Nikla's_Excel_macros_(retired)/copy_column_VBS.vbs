
Dim fpathXlsm
fpathXlsm = wscript.arguments.named.item("fpathXlsm")

Dim fpathGrps
Dim shGrps
Dim colFormatBar
Dim HeaderRow
Dim StartDataRow

fpathGrps = WScript.Arguments.Named.Item("fpathGrps")
shSource = WScript.Arguments.Named.Item("shSource")
shDesti = WScript.Arguments.Named.Item("shDesti")
colToCopy = WScript.Arguments.Named.Item("colToCopy")

HeaderRow = WScript.Arguments.Named.Item("HeaderRow")
StartDataRow = WScript.Arguments.Named.Item("StartDataRow")

Set objExcel = CreateObject("Excel.Application")

Set objWorkbook = objExcel.Workbooks.Open(fpathXlsm)

objExcel.visible = True

objExcel.Application.Run "copy_column_VBA.xlsm!copycol", Cstr(fpathGrps), Cstr(shSource), Cstr(shDesti), Cstr(colToCopy), Clng(HeaderRow), Clng(StartDataRow)

objWorkbook.Close false
objExcel.Quit