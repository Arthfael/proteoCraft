
Dim fpathXlsm
fpathXlsm = wscript.arguments.named.item("fpathXlsm")

Dim fpathGrps
Dim shGrps
Dim colFormatAsNumber
Dim rowStartEnd
Dim countDecimals


fpathGrps = WScript.Arguments.Named.Item("fpathGrps")
shGrps = WScript.Arguments.Named.Item("shGrps")

colFormatAsNumber = WScript.Arguments.Named.Item("colFormatAsNumber")
rowStartEnd = WScript.Arguments.Named.Item("rowStartEnd")
countDecimals = WScript.Arguments.Named.Item("countDecimals")


Set objExcel = CreateObject("Excel.Application")

Set objWorkbook = objExcel.Workbooks.Open(fpathXlsm)

objExcel.visible = True

objExcel.Application.Run "Format_as_num_VBA.xlsm!formatNumeric", Cstr(fpathGrps), Cstr(shGrps), Clng(colFormatAsNumber), Cstr(rowStartEnd), Cbyte(countDecimals)

objWorkbook.Close false
objExcel.Quit