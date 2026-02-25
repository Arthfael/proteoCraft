
Dim fpathXlsm
fpathXlsm = wscript.arguments.named.item("fpathXlsm")

Dim fpathGrps
Dim shGrps
Dim colFormatBar
Dim HeaderRow
Dim StartDataRow

fpathGrps = WScript.Arguments.Named.Item("fpathGrps")
shGrps = WScript.Arguments.Named.Item("shGrps")

colFormatBar = WScript.Arguments.Named.Item("colFormatBar")


StartDataRow = WScript.Arguments.Named.Item("StartDataRow")

Set objExcel = CreateObject("Excel.Application")

Set objWorkbook = objExcel.Workbooks.Open(fpathXlsm)

objExcel.visible = True

objExcel.Application.Run "Format_as_bar_VBA.xlsm!formatBar", Cstr(fpathGrps), Cstr(shGrps), Cstr(colFormatBar), Clng(StartDataRow)

objWorkbook.Close false
objExcel.Quit