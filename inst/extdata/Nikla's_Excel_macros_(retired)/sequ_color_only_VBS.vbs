
Dim fpathXlsm
fpathXlsm = wscript.arguments.named.item("fpathXlsm")

Dim fpathPeptides
Dim shPeptides
Dim colPepID
Dim colPepSequ

Dim fpathProt
Dim shProt
Dim colProtAcc
Dim colProtPepIDs
Dim colProtSequ

Dim fpathGrps
Dim shGrps
Dim colGrpsProtIDs
Dim colGrpsSequ

Dim HeaderRow
Dim StartDataRow

fpathPeptides = WScript.Arguments.Named.Item("fpathPeptides")

colPepID = WScript.Arguments.Named.Item("colPepID")
colPepSequ = WScript.Arguments.Named.Item("colPepSequ")

fpathProt = WScript.Arguments.Named.Item("fpathProt")

colProtAcc = WScript.Arguments.Named.Item("colProtAcc")
colProtPepIDs = WScript.Arguments.Named.Item("colProtPepIDs")
colProtSequ = WScript.Arguments.Named.Item("colProtSequ")

fpathGrps = WScript.Arguments.Named.Item("fpathGrps")
shGrps = WScript.Arguments.Named.Item("shGrps")
colGrpsProtIDs = WScript.Arguments.Named.Item("colGrpsProtIDs")
colGrpsSequ = WScript.Arguments.Named.Item("colGrpsSequ")

HeaderRow = WScript.Arguments.Named.Item("HeaderRow")
StartDataRow = WScript.Arguments.Named.Item("StartDataRow")

Set objExcel = CreateObject("Excel.Application")

Set objWorkbook = objExcel.Workbooks.Open(fpathXlsm)

objExcel.visible = True

objExcel.Application.Run "sequ_color_only_VBA.xlsm!colorSequence", Cstr(fpathPeptides),  Cstr(colPepID), Cstr(colPepSequ), Cstr(fpathProt), Cstr(colProtAcc), Cstr(colProtPepIDs), Cstr(colProtSequ), Cstr(fpathGrps), Cstr(shGrps), Cstr(colGrpsProtIDs), Cstr(colGrpsSequ),  Cstr(HeaderRow), Cstr(StartDataRow)

objWorkbook.Close false
objExcel.Quit