#pragma once
#include <cassert>	
#ifdef PHREEQC_CLASS
#define P_INSTANCE_POINTER1 phreeqc_ptr->
#else
#define P_INSTANCE_POINTER1
#endif
namespace zdg_ui2 {
	using namespace System;
	//using namespace System::ComponentModel;
	using namespace System::Resources;
	using namespace System::Windows::Forms;
	using namespace System::Drawing;
	using namespace System::Threading;
	using namespace ZedGraph;

// Form2 is only used with MULTICHART
	public ref class ChartObj : public System::Object
	{
#ifdef PHREEQC_CLASS
	public: Phreeqc* phreeqc_ptr;
#endif
	public: ChartObject* chartobject_ptr;
	public:	ChartObj(ChartObject* ptr)
			{
				this->chartobject_ptr = ptr;
#ifdef PHREEQC_CLASS
				this->phreeqc_ptr = this->chartobject_ptr->Get_phreeqc();
#endif
			}
	};

	public ref class Form1  : public System::Windows::Forms::Form
	{
	public:	long int tickStart;
	public: Form1 ^myForm;
	public:	Form1()
			{
				InitializeComponent();
				col_use = 0;
				symbol_use = 0;
				Y2 = false;
				phreeqc_done = false;

			}
#ifdef PHREEQC_CLASS
	public:	Form1(ChartObject *ptr)
			{
				this->chartobject_ptr = ptr;
				this->phreeqc_ptr = chartobject_ptr->Get_phreeqc();
				InitializeComponent();
				col_use = 0;
				symbol_use = 0;
				Y2 = false;
				phreeqc_done = false;
				
			}
#else
	public:	Form1(ChartObject *ptr)
			{
				this->chartobject_ptr = ptr;
				InitializeComponent();
				col_use = 0;
				symbol_use = 0;
				Y2 = false;
				phreeqc_done = false;
			}
#endif
			static void ThreadForm(Object^ data)
			{
				ChartObject *ptr = ((ChartObj^)(data))->chartobject_ptr;
				Form1 ^myForm = gcnew Form1(ptr);
				myForm->ShowDialog();
				myForm->~Form1();
			}
	private: bool phreeqc_done;
			 
	private: void SetSize()
			 {
				 zg1->Location = Point( 0, 0 );
				 // Leave a small margin around the outside of the control
				 zg1->Size = System::Drawing::Size( ClientRectangle.Width - 0,
					 ClientRectangle.Height - 0 );
			 }

			System::Void MyFormClosingEventHandler(
					System::Object^ sender, 
				System::Windows::Forms::FormClosingEventArgs ^e)
			{
				ChartObject *chart = this->chartobject_ptr;
				if (chart != NULL) 
				{
					chart->Set_done(true);
					System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
				}
#if defined PHREEQC_CLASS
				this->phreeqc_ptr = NULL;
#endif
				this->chartobject_ptr = NULL;
			}

			 System::Void Form1_Load(System::Object ^sender, System::EventArgs ^e)
			 {
				 CreateGraph( zg1 );
				 SetSize();
			 }

			 System::Void Form1_Resize(System::Object ^sender, System::EventArgs ^e)
			 {
				 SetSize();
			 }

			 static bool LogX, LogY, LogY2;
	private: bool check_neg_log( int i, int i2)
			 {
				 ChartObject *chart = this->chartobject_ptr;
				 if (chart == NULL) return false;
				 std::vector<CurveObject *> &Curves = chart->Get_Curves();
				 if (LogX && chart->Get_axis_scale_x()[4] == 10.0 && 
					 Curves[i]->Get_x()[i2] <= 0)
				 {
					 P_INSTANCE_POINTER1 warning_msg("Obtained x_value <= 0, removing point...");
					 //axis_scale_x[4] = NA; /* if reverting to linear... */
					 //LogX = false;
					 return true;
				 }
				 if (Curves[i]->Get_y()[i2] <= 0 && 
					 (chart->Get_axis_scale_y()[4] == 10.0 || 
					 chart->Get_axis_scale_y2()[4] == 10.0))
				 {
					 if (Curves[i]->Get_y_axis() == 2 && LogY2)
					 {
						 P_INSTANCE_POINTER1 warning_msg("Obtained sy_value <= 0, removing point......");
						 //axis_scale_y2[4] = NA;
						 //LogY2 = false;
						 return true;
					 }
					 else if (LogY)
					 {
						 P_INSTANCE_POINTER1 warning_msg("Obtained y_value <= 0, removing point......");
						 //axis_scale_y[4] = NA;
						 //LogY = false;
						 return true;
					 }
				 }
				 return false;
			 }

	private: PointPairList ^list;
			 int col_use, symbol_use;
			 bool Y2;
			 static cli::array<String^> ^ColorList = {"Red", "Green", "Blue", "Orange", "Magenta", "Yellow", "Black" };

			 void DefineCurves(GraphPane ^myPane, int init)
			 {
				 ChartObject *chart = this->chartobject_ptr;
				 if (chart == NULL) 
				 {
					 return;
				 }
 				 std::vector<CurveObject *> Curves; 
				 size_t i;
				 for (i = 0; i < chart->Get_CurvesCSV().size(); i++)
				 {
					 Curves.push_back(chart->Get_CurvesCSV()[i]);
				 }
				 for (i = 0; i < chart->Get_Curves().size(); i++)
				 {
					 Curves.push_back(chart->Get_Curves()[i]);
				 }

				 chart->Set_curve_added(false);

				 // Set the titles and axis labels
				 myPane->Title->Text = gcnew String(chart->Get_chart_title().c_str());
				 if (chart->Get_axis_titles().size() > 0)
					 myPane->XAxis->Title->Text = gcnew String(chart->Get_axis_titles()[0].c_str());
				 if (chart->Get_axis_titles().size() > 1)
					 myPane->YAxis->Title->Text = gcnew String(chart->Get_axis_titles()[1].c_str());
				 if (chart->Get_axis_titles().size() > 2)
					 myPane->Y2Axis->Title->Text = gcnew String(chart->Get_axis_titles()[2].c_str());

				 LineItem ^myCurve;

				 Color col;

				 String ^s_t;
				 if (chart->Get_axis_scale_x()[4] == 10.0) LogX = true;
				 else LogX = false;
				 if (chart->Get_axis_scale_y()[4] == 10.0) LogY = true;
				 else LogY = false;
				 if (chart->Get_axis_scale_y2()[4] == 10.0) LogY2 = true;
				 else LogY2 = false;

				 //Rewrite all curves
				 zg1->GraphPane->CurveList->Clear();
				 for (size_t i = 0; i < Curves.size(); i++)
				 {
					 // even curves with no data
					 //if (Curves[i]->Get_x().size() == 0) continue;
					 list = gcnew PointPairList();
					 if (Curves[i]->Get_y_axis() == 2)
						 Y2 = true;
					 else
						 Y2 = false;
					 for (int i2 = 0; (i2 < (int) Curves[i]->Get_x().size()); i2++)
					 {
						 if ((LogX && Curves[i]->Get_x()[i2] <=0)
							 || (LogY && !Y2 && Curves[i]->Get_y()[i2] <=0)
							 || (LogY2 && Y2 && Curves[i]->Get_y()[i2] <=0))
							 continue;
						 else
							 list->Add( Curves[i]->Get_x()[i2], Curves[i]->Get_y()[i2] );
					 }

					 col = Color::FromName(gcnew String(Curves[i]->Get_color().c_str()));
					 if (!col.IsKnownColor)
					 {
						 col = Color::FromName(ColorList[col_use]);
						 std::string newcol;
						 ToString(col.ToString(), newcol);
						 Utilities::replace("Color [","",newcol);
						 Utilities::replace("]","",newcol);
						 Curves[i]->Set_color(newcol);

					 }
					 if (++col_use > 6) col_use = 0;

					 SymbolType symb = chart->Return_SymbolType
						 (Curves[i]->Get_symbol());

					 // id
					 s_t = gcnew String(Curves[i]->Get_id().c_str());

					 // Curve with no points is invisible
					 if (Curves[i]->Get_x().size() == 0) 
					 {
						 s_t = "";
					 }

					 // Add curve to chart
					 myCurve = myPane->AddCurve( s_t, list, col, symb );

					 // Curve with no points is invisible
					 if (Curves[i]->Get_x().size() == 0) 
					 {
						 myCurve->IsVisible = false;
					 }

					 if (Curves[i]->Get_line_w() > 0.0)
						 myCurve->Line->Width = (float) Curves[i]->Get_line_w();
					 else
						 myCurve->Line->IsVisible = false;
					 /* hmm... dash/dot don't display well */
					 //myCurve->Line->Style = System::Drawing::Drawing2D::DashStyle::Dot;
					 myCurve->Symbol->Fill = gcnew Fill( Color::FromName("White") );
					 if (Curves[i]->Get_symbol_size() > 0.0)
						 myCurve->Symbol->Size = (float) Curves[i]->Get_symbol_size();
					 else
						 myCurve->Symbol->IsVisible = false;
					 myCurve->Symbol->Border->Width = (float) Curves[i]->Get_line_w();
					 if (Y2)
						 myCurve->IsY2Axis = true;
					 delete list;
				 }

				 if (Y2)
					 myPane->Legend->Position = ZedGraph::LegendPos::TopCenter;
				 else
					 myPane->Legend->Position = ZedGraph::LegendPos::Right;
				 myPane->Legend->FontSpec->Size = 12;
				 myPane->Legend->FontSpec->IsBold = false;

				 // Show the x axis grid
				 myPane->XAxis->MajorGrid->IsVisible = true;
				 if (fabs(chart->Get_axis_scale_x()[0] - NA) > 1e-3)
					 myPane->XAxis->Scale->Min = chart->Get_axis_scale_x()[0];
				 else
					 myPane->XAxis->Scale->MinAuto = true;
				 if (fabs(chart->Get_axis_scale_x()[1] - NA) > 1e-3)
					 myPane->XAxis->Scale->Max = chart->Get_axis_scale_x()[1];
				 else
					 myPane->XAxis->Scale->MaxAuto = true;
				 if (fabs(chart->Get_axis_scale_x()[2] - NA) > 1e-3)
					 myPane->XAxis->Scale->MajorStep = chart->Get_axis_scale_x()[2];
				 else
					 myPane->XAxis->Scale->MajorStepAuto = true;
				 if (fabs(chart->Get_axis_scale_x()[3] - NA) > 1e-3)
				 {
					 myPane->XAxis->Scale->MinorStep = chart->Get_axis_scale_x()[3];
					 if (chart->Get_axis_scale_x()[3] == 0.0)
						 // remove minor tics
						 myPane->XAxis->MinorTic->Size = 0;
				 }
				 else
					 myPane->XAxis->Scale->MinorStepAuto = true;
				 if (chart->Get_axis_scale_x()[4] == 10.0)
					 myPane->XAxis->Type = AxisType::Log;

				 // Make the Y axis scale red
				 // myPane->YAxis->Scale->FontSpec->FontColor = Color::Red;
				 // myPane->YAxis->Title->FontSpec->FontColor = Color::Red;
				 // turn off the opposite tics so the Y tics don't show up on the Y2 axis
				 if (Y2)
				 {
					 myPane->YAxis->MajorTic->IsOpposite = false;
					 myPane->YAxis->MinorTic->IsOpposite = false;
				 }
				 // Don't display the Y zero line
				 myPane->YAxis->MajorGrid->IsZeroLine = false;
				 // Align the Y axis labels so they are flush to the axis
				 myPane->YAxis->Scale->Align = AlignP::Inside;
				 myPane->YAxis->MajorGrid->IsVisible = true;
				 if (fabs(chart->Get_axis_scale_y()[0] - NA) > 1e-3)
					 myPane->YAxis->Scale->Min = chart->Get_axis_scale_y()[0];
				 else
					 myPane->YAxis->Scale->MinAuto = true;
				 if (fabs(chart->Get_axis_scale_y()[1] - NA) > 1e-3)
					 myPane->YAxis->Scale->Max = chart->Get_axis_scale_y()[1];
				 else
					 myPane->YAxis->Scale->MaxAuto = true;
				 if (fabs(chart->Get_axis_scale_y()[2] - NA) > 1e-3)
					 myPane->YAxis->Scale->MajorStep = chart->Get_axis_scale_y()[2];
				 else
					 myPane->YAxis->Scale->MajorStepAuto = true;
				 if (fabs(chart->Get_axis_scale_y()[3] - NA) > 1e-3)
				 {
					 myPane->YAxis->Scale->MinorStep = chart->Get_axis_scale_y()[3];
					 if (chart->Get_axis_scale_y()[3] == 0.0)
						 // remove minor tics
						 myPane->YAxis->MinorTic->Size = 0;
				 }
				 else
					 myPane->YAxis->Scale->MinorStepAuto = true;
				 if (chart->Get_axis_scale_y()[4] == 10.0)
					 myPane->YAxis->Type = AxisType::Log;

				 // Enable the Y2 axis display
				 if (Y2)
				 {
					 myPane->Y2Axis->IsVisible = true;
					 // Make the Y2 axis scale blue
					 // myPane->Y2Axis->Scale->FontSpec->FontColor = Color::Blue;
					 // myPane->Y2Axis->Title->FontSpec->FontColor = Color::Blue;
					 // turn off the opposite tics so the Y2 tics don't show up on the Y axis
					 myPane->Y2Axis->MajorTic->IsOpposite = false;
					 myPane->Y2Axis->MinorTic->IsOpposite = false;
					 // Don't display the Y2 axis grid lines
					 myPane->Y2Axis->MajorGrid->IsVisible = false;
					 // Align the Y2 axis labels so they are flush to the axis
					 myPane->Y2Axis->Scale->Align = AlignP::Inside;

					 if (fabs(chart->Get_axis_scale_y2()[0] - NA) > 1e-3)
						 myPane->Y2Axis->Scale->Min = chart->Get_axis_scale_y2()[0];
					 else
						 myPane->Y2Axis->Scale->MinAuto = true;
					 if (fabs(chart->Get_axis_scale_y2()[1] - NA) > 1e-3)
						 myPane->Y2Axis->Scale->Max = chart->Get_axis_scale_y2()[1];
					 else
						 myPane->Y2Axis->Scale->MaxAuto = true;
					 if (fabs(chart->Get_axis_scale_y2()[2] - NA) > 1e-3)
						 myPane->Y2Axis->Scale->MajorStep = chart->Get_axis_scale_y2()[2];
					 else
						 myPane->Y2Axis->Scale->MajorStepAuto = true;
					 if (fabs(chart->Get_axis_scale_y2()[3] - NA) > 1e-3)
					 {
						 myPane->Y2Axis->Scale->MinorStep = chart->Get_axis_scale_y2()[3];
						 if (chart->Get_axis_scale_y2()[3] == 0.0)
							 // remove minor tics
							 myPane->Y2Axis->MinorTic->Size = 0;
					 }
					 else
						 myPane->Y2Axis->Scale->MinorStepAuto = true;
					 if (chart->Get_axis_scale_y2()[4] == 10.0)
						 myPane->Y2Axis->Type = AxisType::Log;
				 }

				 // Fill the axis background with a gradient
				 //myPane->Chart->Fill = gcnew Fill( Color::White, Color::LightYellow, 45.0f ); /* FromArgb(255, 255, 224) */
				 myPane->Chart->Fill = gcnew Fill( Color::White, Color::FromArgb(255, 255, 230), 45.0f );

				 // Make sure auto scale, Refresh
				 zg1->AxisChange();
				 zg1->Refresh();

			 }

	public: void CreateGraph( ZedGraphControl ^z1 )	{
				// Get a reference to the GraphPane instance in the ZedGraphControl
				GraphPane ^myPane = z1->GraphPane;

				// lock thread
				while (0 != System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 1)) 
					System::Threading::Thread::Sleep(1);

				DefineCurves(myPane, 0);

				// Add text boxes with instructions...
				TextObj ^text;
				text = gcnew TextObj(
					L" Click right mouse for options...",
					0.01f, 0.99f, CoordType::PaneFraction, AlignH::Left, AlignV::Bottom );
				text->FontSpec->StringAlignment = StringAlignment::Near;
				text->FontSpec->Size = 10;
				text->FontSpec->FontColor = Color::Red;
				myPane->GraphObjList->Add( text );
				text = gcnew TextObj(
					L" Press Alt + F4 to quit",
					0.81f, 0.99f, CoordType::PaneFraction, AlignH::Left, AlignV::Bottom );
				text->FontSpec->StringAlignment = StringAlignment::Near;
				text->FontSpec->Size = 10;
				text->FontSpec->FontColor = Color::Red;
				myPane->GraphObjList->Add( text );

				// Enable scrollbars if needed...
				/*z1->IsShowHScrollBar = true;
				z1->IsShowVScrollBar = true;
				z1->IsAutoScrollRange = true;
				z1->IsScrollY2 = true;*/

				// OPTIONAL: Show tooltips when the mouse hovers over a point
				z1->IsShowPointValues = false;
				z1->PointValueEvent += gcnew ZedGraphControl::PointValueHandler( this,
					&Form1::MyPointValueHandler );

				// OPTIONAL: Add a custom context menu item
				z1->ContextMenuBuilder += gcnew	ZedGraphControl::ContextMenuBuilderEventHandler(
					this, &Form1::MyContextMenuBuilder );

				// OPTIONAL: Handle the Zoom Event
				z1->ZoomEvent += gcnew ZedGraphControl::ZoomEventHandler( this,
					&Form1::MyZoomEvent );

				// Size the control to fit the window
				SetSize();

				// Tell ZedGraph to calculate the axis ranges
				// Note that you MUST call this after enabling IsAutoScrollRange, since AxisChange() sets
				// up the proper scrolling parameters

				z1->AxisChange();
				// Make sure the Graph gets redrawn
				z1->Invalidate();
				timer1->Interval = this->chartobject_ptr->Get_update_time_chart();
				timer1->Enabled = true;
				timer1->Start();

				tickStart = Environment::TickCount;

				//unlock thread
				System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
			}

			/// <summary>
			/// Display customized tooltips when the mouse hovers over a point
			/// </summary>
			System::String ^MyPointValueHandler( ZedGraphControl ^control, GraphPane ^pane,
				CurveItem ^curve, int iPt ) {
					// Get the PointPair that is under the mouse
					PointPair pt = curve[iPt];
					return curve->Label->Text + " is " + pt.Y.ToString( "e3" ) + " units at X = " + pt.X.ToString( "e3" );
			}

			// Add some explanation to the menu..
			void MyContextMenuBuilder( ZedGraphControl ^control,
				System::Windows::Forms::ContextMenuStrip ^menuStrip,
				Point mousePt,
				ZedGraphControl::ContextMenuObjectState objState ) {
					ToolStripMenuItem ^item = gcnew ToolStripMenuItem();
					item->Text = L"Zoom: left mouse + drag\nPan: middle mouse + drag";
					menuStrip->Items->Insert(5, item );

					menuStrip->Items->RemoveAt(0);
					ToolStripMenuItem ^item2 = gcnew ToolStripMenuItem();
					item2->Text = L"Save Data to File \'curves.u_g\'";
					item2->Click += gcnew System::EventHandler(this, &zdg_ui2::Form1::SaveCurves );
					menuStrip->Items->Insert(0, item2 );

			}

			void form_error_msg( std::string estring )
			{
				if (this->chartobject_ptr != NULL)
				{
					P_INSTANCE_POINTER1 error_msg(estring.c_str(), CONTINUE);
				}
				else
				{
					std::cerr << "ERROR: " << estring << std::endl;
				}
			}

			void SaveCurves( System::Object ^sender, System::EventArgs ^e )
			{
				std::string file_name = "curves.u_g";
						 // Get the graph curves...
				std::ofstream f_out(file_name.c_str(), std::ifstream::out);

				if (!f_out.is_open())
				{
					std::ostringstream estream;
					estream << "Could not open csv file for USER_GRAPH " << file_name;
					form_error_msg(estream.str());
					return;
				}

				// write headings
				size_t max_points = 0; 
				f_out.precision(4);
				for (int i = 0; i < zg1->GraphPane->CurveList->Count; i++) 
				{
					LineItem  ^curve;
					curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
					// Get the PointPairList
					IPointListEdit  ^ip = (IPointListEdit^) curve->Points;

					// Calculate max_points
					if ((size_t) ip->Count > max_points)
						max_points = ip->Count;

					// write headers
					std::string s_std;
					ToString(curve->Label->Text, s_std);
					f_out.width(12);
					f_out << "x" << "\t";
					f_out.width(12);
					if (s_std.size() > 0) 
					{
						f_out << s_std << "\t";
					}
					else
					{
						f_out << "y" << "\t";
					}
				}

				f_out << std::endl;

				// write data
				size_t i2 = 0;
				f_out << std::scientific;
				f_out.precision(4);

				while (i2 < max_points)
				{
					for (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
					{
						LineItem  ^curve;
						curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
						// Get the PointPairList
						IPointListEdit  ^ip = (IPointListEdit^) curve->Points;					
						if (i2 < (size_t) ip->Count)
						{
							//double x = ip[i]->X;
							f_out.width(12);
							f_out << ip[i2]->X << "\t";
							f_out.width(12);
							f_out << ip[i2]->Y << "\t";
						}
						else if (i2 < max_points)
						{
							f_out.width(13);
							f_out << "\t";
							f_out.width(13);
							f_out << "\t";
						}
					}
					f_out << std::endl;
					i2++;
				}
				f_out.close();
				return;	
			}


			// Respond to a Zoom Event
			void MyZoomEvent( ZedGraphControl ^control, ZoomState ^oldState, ZoomState ^newState )
			{
				// Here we get notification everytime the user zooms
			}
   private: void timer1_Tick(System::Object ^sender, System::EventArgs ^e )
			{
				LineItem  ^curve;
				ChartObject *chart = this->chartobject_ptr;
				if (chart == NULL) return;

				//lock for thread
				while (0 != System::Threading::Interlocked::Exchange(chart->usingResource, 1)) 
					System::Threading::Thread::Sleep(1);

				if (this->chartobject_ptr->Get_curve_added())
				{
					DefineCurves(zg1->GraphPane, zg1->GraphPane->CurveList->Count);
				}
				else if (this->chartobject_ptr->Get_point_added())
				{

					// Make list of curves
					std::vector<CurveObject *> Curves; 
					size_t j;
					for (j = 0; j < chart->Get_CurvesCSV().size(); j++)
					{
						Curves.push_back(chart->Get_CurvesCSV()[j]);
					}
					for (j = 0; j < chart->Get_Curves().size(); j++)
					{
						Curves.push_back(chart->Get_Curves()[j]);
					}
					// Add points to curves ...
					for (int i = 0; i < zg1->GraphPane->CurveList->Count; i++) 
					{
						curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
						// Get the PointPairList
						IPointListEdit  ^ip = (IPointListEdit^) curve->Points;
						if ((size_t) ip->Count < Curves[i]->Get_x().size())
						{
							if (Curves[i]->Get_y_axis() == 2)
								Y2 = true;
							else
								Y2 = false;
							for ( size_t i2 = ip->Count; i2 < Curves[i]->Get_x().size(); i2++ )
							{
								if ((LogX && Curves[i]->Get_x()[i2] <=0)
									|| (LogY && !Y2 && Curves[i]->Get_y()[i2] <=0)
									|| (LogY2 && Y2 && Curves[i]->Get_y()[i2] <=0))
									continue;
								else
									ip->Add(Curves[i]->Get_x()[i2], Curves[i]->Get_y()[i2] );
							}
						}
					}
					// Add points to curves ...

					//size_t i, k;
					//k = 0;
					//for (i = 0; i < Curves.size(); i++) 
					//{
					//	if (Curves[i]->Get_x().size() == 0) continue;
					//	curve =  (LineItem ^) zg1->GraphPane->CurveList[k++];
					//	// Get the PointPairList
					//	IPointListEdit  ^ip = (IPointListEdit^) curve->Points;
					//	if ((size_t) ip->Count < Curves[i]->Get_x().size())
					//	{
					//		if (Curves[i]->Get_y_axis() == 2)
					//			Y2 = true;
					//		else
					//			Y2 = false;
					//		for ( size_t i2 = ip->Count; i2 < Curves[i]->Get_x().size(); i2++ )
					//		{
					//			if ((LogX && Curves[i]->Get_x()[i2] <=0)
					//				|| (LogY && !Y2 && Curves[i]->Get_y()[i2] <=0)
					//				|| (LogY2 && Y2 && Curves[i]->Get_y()[i2] <=0))
					//				continue;
					//			else
					//				ip->Add(Curves[i]->Get_x()[i2], Curves[i]->Get_y()[i2] );
					//		}
					//	}
					//}
					/* explicitly reset the max in case of log scale, zedgraphs doesn't do this... */
					if ((fabs(chart->Get_axis_scale_x()[1] - NA) < 1e-3) && zg1->GraphPane->XAxis->Type == AxisType::Log)
					{
						double max = -1e99;
						for  (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
						{
							if (Curves[i]->Get_x()[Curves[i]->Get_x().size() - 1] > max)
								max = Curves[i]->Get_x()[Curves[i]->Get_x().size() - 1];
						}
						max += pow(10.0, log10(max / 3));
						zg1->GraphPane->XAxis->Scale->Max = max;
					}
					if ((fabs(chart->Get_axis_scale_y()[1] - NA) < 1e-3) && zg1->GraphPane->YAxis->Type == AxisType::Log)
					{
						double max = -1e99;
						for  (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
						{
							curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
							if (curve->IsY2Axis) continue;
							if (Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1] > max)
								max = Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1];
						}
						max += pow(10.0, log10(max / 3));
						zg1->GraphPane->YAxis->Scale->Max = max;
					}
					if ((fabs(chart->Get_axis_scale_y2()[1] - NA) < 1e-3) && zg1->GraphPane->Y2Axis->Type == AxisType::Log)
					{
						double max = -1e99;
						for  (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
						{
							curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
							if (!curve->IsY2Axis) continue;
							if (Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1] > max)
								max = Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1];
						}
						max += pow(10.0, log10(max / 3));
						zg1->GraphPane->Y2Axis->Scale->Max = max;
					}
					zg1->AxisChange();
					zg1->Refresh();
				}

				chart->Set_point_added(false);
				if (chart->Get_end_timer())
				{
					timer1->Stop();
					chart->Set_done(true);
					phreeqc_done = true;

					{
						//// Debugging check
						//std::vector<CurveObject *> Curves = chart->Get_CurvesCSV(); 
						//size_t i;
						//for (i = 0; i < chart->Get_Curves().size(); i++)
						//{
						//	Curves.push_back(chart->Get_Curves()[i]);
						//}
						//for (i = 0; i < (size_t) zg1->GraphPane->CurveList->Count; i++) 
						//{
						//	if (zg1->GraphPane->CurveList[i]->Points->Count !=
						//		Curves[i]->Get_x().size())
						//	{
						//		std::cerr << "Form: " << i << "\t" << zg1->GraphPane->CurveList[i]->Points->Count << std::endl;
						//		std::cerr << "Curves: " << i << "\t" << Curves[i]->Get_x().size() << std::endl;
						//		//form_error_msg("Did not plot all points. Why?");
						//	}
						//	assert (zg1->GraphPane->CurveList[i]->Points->Count ==
						//		Curves[i]->Get_x().size());
						//}
					}
					//unlock thread before setting chartobject_ptr to NULL
					System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
#if defined PHREEQC_CLASS
					this->phreeqc_ptr = NULL;
#endif
					this->chartobject_ptr = NULL;
					return;
				}

				//unlock thread
				System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
				//tickStart = Environment::TickCount;
				return;
			}

	private: void timer1_Tick_old(System::Object ^sender, System::EventArgs ^e )
			 {
				 LineItem  ^curve;
				 ChartObject *chart = this->chartobject_ptr;
				 if (chart == NULL) return;
				 //std::cerr << "timer1_Tick." << std::endl;
				 //lock for thread
				 while (0 != System::Threading::Interlocked::Exchange(chart->usingResource, 1)) 
					 System::Threading::Thread::Sleep(1);
				 
				 std::vector<CurveObject *> Curves; 
				 size_t i;
				 for (i = 0; i < chart->Get_CurvesCSV().size(); i++)
				 {
					 Curves.push_back(chart->Get_CurvesCSV()[i]);
				 }
				 for (i = 0; i < chart->Get_Curves().size(); i++)
				 {
					 Curves.push_back(chart->Get_Curves()[i]);
				 }

				 if ( ((Environment::TickCount - tickStart ) > this->chartobject_ptr->Get_update_time_chart())
					 || chart->Get_end_timer() ) {
					if (this->chartobject_ptr->Get_curve_added())
					 {
						 DefineCurves(zg1->GraphPane, zg1->GraphPane->CurveList->Count);

					 }
					 else if (this->chartobject_ptr->Get_point_added())
					 {
						 // Add points to curves ...
						 for (int i = 0; i < zg1->GraphPane->CurveList->Count; i++) 
						 {
							 curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
							 // Get the PointPairList
							 IPointListEdit  ^ip = (IPointListEdit^) curve->Points;
							 if ((size_t) ip->Count < Curves[i]->Get_x().size())
							 {

								 if (Curves[i]->Get_y_axis() == 2)
									 Y2 = true;
								 else
									 Y2 = false;
								 for ( size_t i2 = ip->Count; i2 < Curves[i]->Get_x().size(); i2++ )
								 {
									 if ((LogX && Curves[i]->Get_x()[i2] <=0)
										 || (LogY && !Y2 && Curves[i]->Get_y()[i2] <=0)
										 || (LogY2 && Y2 && Curves[i]->Get_y()[i2] <=0))
										 continue;
									 else
										 ip->Add(Curves[i]->Get_x()[i2], Curves[i]->Get_y()[i2] );
								 }
							 }
						 }

						 /* explicitly reset the max in case of log scale, zedgraphs doesn't do this... */
						 if ((fabs(chart->Get_axis_scale_x()[1] - NA) < 1e-3) && zg1->GraphPane->XAxis->Type == AxisType::Log)
						 {
							 double max = -1e99;
							 for  (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
							 {
								 if (Curves[i]->Get_x()[Curves[i]->Get_x().size() - 1] > max)
									 max = Curves[i]->Get_x()[Curves[i]->Get_x().size() - 1];
							 }
							 max += pow(10.0, log10(max / 3));
							 zg1->GraphPane->XAxis->Scale->Max = max;
						 }
						 if ((fabs(chart->Get_axis_scale_y()[1] - NA) < 1e-3) && zg1->GraphPane->YAxis->Type == AxisType::Log)
						 {
							 double max = -1e99;
							 for  (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
							 {
								 curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
								 if (curve->IsY2Axis) continue;
								 if (Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1] > max)
									 max = Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1];
							 }
							 max += pow(10.0, log10(max / 3));
							 zg1->GraphPane->YAxis->Scale->Max = max;
						 }
						 if ((fabs(chart->Get_axis_scale_y2()[1] - NA) < 1e-3) && zg1->GraphPane->Y2Axis->Type == AxisType::Log)
						 {
							 double max = -1e99;
							 for  (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
							 {
								 curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
								 if (!curve->IsY2Axis) continue;
								 if (Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1] > max)
									 max = Curves[i]->Get_y()[Curves[i]->Get_y().size() - 1];
							 }
							 max += pow(10.0, log10(max / 3));
							 zg1->GraphPane->Y2Axis->Scale->Max = max;
						 }

						 zg1->AxisChange();
						 zg1->Refresh();
						 tickStart = Environment::TickCount;
					 }
				 }
				 chart->Set_point_added(false);
				 if (chart->Get_end_timer())
				 {
					 //std::cerr << "Form got end_timer message." << std::endl;
					 zg1->Refresh();
					 timer1->Stop();
					 chart->Set_done(true);
					 phreeqc_done = true;
					 //unlock thread
					 System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
#if defined PHREEQC_CLASS
					this->phreeqc_ptr = NULL;
#endif
					this->chartobject_ptr = NULL;
					//std::cerr << "Form released thread, pointers null." << std::endl;
				 }
				 else
				 {
					 //unlock thread
					 System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
				 }
				 return;
			 }
			 void ToString(System::String^ src, std::string& dest)
			 {
				 using namespace System::Runtime::InteropServices;
				 const char* chars = (const char*)(Marshal::StringToHGlobalAnsi(src)).ToPointer();
				 dest = chars;
				 Marshal::FreeHGlobal(IntPtr((void*)chars));
			 }
			 ~Form1() {
				 if (this->zg1) delete zg1;
				 //if (this->timer1) delete timer1);
				 if (components) {
					 delete components;
				 }
			 }
	public: ZedGraph::ZedGraphControl ^zg1;
	private: System::Windows::Forms::Timer ^timer1;
	private: System::ComponentModel::Container ^components;
#ifdef PHREEQC_CLASS
	private: Phreeqc * phreeqc_ptr;
#endif
			 ChartObject * chartobject_ptr;

	public:
#pragma region Windows Form Designer generated code
		/// <summary>
		/// Required method for Designer support - do not modify
		/// the contents of this method with the code editor.
		/// </summary>
		void InitializeComponent()
		{
			this->components = (gcnew System::ComponentModel::Container());
			this->zg1 = (gcnew ZedGraph::ZedGraphControl());
			this->timer1 = (gcnew System::Windows::Forms::Timer( this->components ));
			this->SuspendLayout();
			// 
			// zg1
			// 
			this->zg1->Location = System::Drawing::Point(12, 12);
			this->zg1->Name = L"zg1";
			this->zg1->ScrollGrace = 0;
			this->zg1->ScrollMaxX = 0;
			this->zg1->ScrollMaxY = 0;
			this->zg1->ScrollMaxY2 = 0;
			this->zg1->ScrollMinX = 0;
			this->zg1->ScrollMinY = 0;
			this->zg1->ScrollMinY2 = 0;
			this->zg1->Size = System::Drawing::Size(this->chartobject_ptr->Get_PanelWidth() - 2 * 12, chartobject_ptr->Get_PanelHeight() - 2 * 12);
			this->zg1->TabIndex = 0;
			this->timer1->Tick += gcnew System::EventHandler( this, &Form1::timer1_Tick );
			// 
			// Form1
			// 
			this->AutoScaleDimensions = System::Drawing::SizeF(6, 13);
			this->AutoScaleMode = System::Windows::Forms::AutoScaleMode::Font;
			this->AutoValidate = System::Windows::Forms::AutoValidate::EnablePreventFocusChange;
			this->ClientSize = System::Drawing::Size(this->chartobject_ptr->Get_PanelWidth(), chartobject_ptr->Get_PanelHeight());
			this->Controls->Add(this->zg1);
			this->Name = L"Form1";
			this->StartPosition = System::Windows::Forms::FormStartPosition::WindowsDefaultLocation;//:CenterScreen;
			this->Text = L"PHREEQC chart";
			this->TopMost = false;
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(Form1::typeid));
			try
			{
				this->Icon = (cli::safe_cast<System::Drawing::Icon^  >(resources->GetObject(L"$this.Icon")));
			}
			catch (...)
			{
			}

			this->FormClosing += gcnew System::Windows::Forms::FormClosingEventHandler(this, &Form1::MyFormClosingEventHandler);
			this->Load += gcnew System::EventHandler(this, &Form1::Form1_Load);
			this->Resize += gcnew System::EventHandler(this, &Form1::Form1_Resize);
			this->ResumeLayout(false);
		}
#pragma endregion


	};
}