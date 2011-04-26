#pragma once
#include <cassert>	
#ifdef PHREEQC_CLASS
//#define P_INSTANCE p_instance
//#define P_INSTANCE_COMMA p_instance,
#define P_INSTANCE_POINTER phreeqc_ptr->
//#define PHREEQC_PTR_ARG Phreeqc *p_instance
//#define PHREEQC_PTR_ARG_COMMA Phreeqc *p_instance,
#else
#define P_INSTANCE_POINTER
#endif
namespace zdg_ui2 {
	using namespace System;
	//using namespace System::ComponentModel;
	using namespace System::Resources;
	using namespace System::Windows::Forms;
	using namespace System::Drawing;
	using namespace System::Threading;
	using namespace ZedGraph;

#ifdef PHREEQC_CLASS
	public ref class PhreeqcObj : public System::Object
	{
	public: Phreeqc* phreeqc_ptr;
	public:	PhreeqcObj(Phreeqc* ptr)
	{
		this->phreeqc_ptr = ptr;
	}
	};
#endif
#ifdef MULTICHART
	public ref class ChartObj : public System::Object
	{
	public: ChartObject* chartobject_ptr;
	public:	ChartObj(ChartObject* ptr)
	{
		this->chartobject_ptr = ptr;
	}
	};
#endif

	public ref class Form1  : public System::Windows::Forms::Form
	{
	public:	long int tickStart;
	public: Form1 ^myForm;
	public:	Form1()
	{
#ifdef PHREEQC_CLASS
		this->phreeqc_ptr = NULL;
#endif
#ifdef MULTICHART
		this->chartobject_ptr = NULL;
#endif
		InitializeComponent();
	}
#ifdef PHREEQC_CLASS
	public:	Form1(Phreeqc *ptr)
	{
		this->phreeqc_ptr = ptr;
		InitializeComponent();
	}
#endif
	public:	Form1(ChartObject *ptr)
	{
		this->chartobject_ptr = ptr;
		InitializeComponent();
	}
#ifdef PHREEQC_CLASS
	static void ThreadForm(Object^ data)
	{
		Phreeqc *ptr = ((PhreeqcObj^)(data))->phreeqc_ptr;
		ptr->u_g_active = true;
		Form1 ^myForm = gcnew Form1(ptr);
		myForm->ShowDialog();
		ptr->u_g_active = false;
	}
#elif defined(MULTICHART)
	static void ThreadForm(Object^ data)
	{

		ChartObject *ptr = ((ChartObj^)(data))->chartobject_ptr;
		Form1 ^myForm = gcnew Form1(ptr);
		myForm->ShowDialog();
		myForm->~Form1();

	}
#else
	static void ThreadForm()
	{
		Form1 ^myForm = gcnew Form1();
		myForm->ShowDialog();

	}
#endif
	//void ThreadThis(Object^ data)
	//{
	//	this->ShowDialog();
	//	assert(false);
	//}
	private: void SetSize()
	{
		zg1->Location = Point( 0, 0 );
			// Leave a small margin around the outside of the control
		zg1->Size = System::Drawing::Size( ClientRectangle.Width - 0,
				ClientRectangle.Height - 0 );
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
			std::vector<CurveObject> &Curves = chart->Get_Curves();
			if (LogX && chart->Get_axis_scale_x()[4] == 10.0 && 
				Curves[i].Get_x()[i2] <= 0)
			{
				P_INSTANCE_POINTER warning_msg("Obtained x_value <= 0, removing point...");
				//axis_scale_x[4] = NA; /* if reverting to linear... */
				//LogX = false;
				return true;
			}
			if (Curves[i].Get_y()[i2] <= 0 && 
			   (chart->Get_axis_scale_y()[4] == 10.0 || 
			   chart->Get_axis_scale_y2()[4] == 10.0))
			{
				if (Curves[i].Get_y_axis() == 2 && LogY2)
				{
					P_INSTANCE_POINTER warning_msg("Obtained sy_value <= 0, removing point......");
					//axis_scale_y2[4] = NA;
					//LogY2 = false;
				return true;
				}
				else if (LogY)
				{
					P_INSTANCE_POINTER warning_msg("Obtained y_value <= 0, removing point......");
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
			if (chart == NULL) return;

			std::vector<CurveObject> Curves = chart->Get_CurvesCSV();
			size_t i;
			for (i = 0; i < chart->Get_CurvesPrevious().size(); i++)
			{
				Curves.push_back(chart->Get_CurvesPrevious()[i]);
			}
			for (i = 0; i < chart->Get_Curves().size(); i++)
			{
				Curves.push_back(chart->Get_Curves()[i]);
			}

			chart->Get_ncurves_changed()[0] = 0;

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
				if (Curves[i].Get_x().size() == 0) continue;
				list = gcnew PointPairList();
				if (Curves[i].Get_y_axis() == 2)
					Y2 = true;
				else
					Y2 = false;
				for (int i2 = 0; i2 < (int) Curves[i].Get_x().size(); i2++)
				{
					if ((LogX && Curves[i].Get_x()[i2] <=0)
						|| (LogY && !Y2 && Curves[i].Get_y()[i2] <=0)
						|| (LogY2 && Y2 && Curves[i].Get_y()[i2] <=0))
						continue;
					else
						list->Add( Curves[i].Get_x()[i2], 
						Curves[i].Get_y()[i2] );
				}

				// Get legal color
				if (strlen(Curves[i].Get_color().c_str()) > 0) {
					col = Color::FromName(gcnew String(Curves[i].Get_color().c_str()));
					if (!col.IsKnownColor)
					{
						col = Color::FromName(ColorList[col_use]);
						std::string newcol;
						ToString(col.ToString(), newcol);
						Curves[i].Set_color(newcol);
					}
				}
				else 
				{
					col = Color::FromName(ColorList[col_use]);
					std::string newcol;
					ToString(col.ToString(), newcol);
					Curves[i].Set_color(newcol);
				}
				if (++col_use > 6) col_use = 0;


				SymbolType symb = chart->Return_SymbolType
					(Curves[i].Get_symbol());

				s_t = gcnew String(Curves[i].Get_id().c_str());

				// Add curve to chart
				myCurve = myPane->AddCurve( s_t, list, col, symb );

				if (Curves[i].Get_line_w() > 0.0)
					myCurve->Line->Width = (float) Curves[i].Get_line_w();
				else
					myCurve->Line->IsVisible = false;
				/* hmm... dash/dot don't display well */
				//myCurve->Line->Style = System::Drawing::Drawing2D::DashStyle::Dot;
				myCurve->Symbol->Fill = gcnew Fill( Color::FromName("White") );
				if (Curves[i].Get_symbol_size() > 0.0)
					myCurve->Symbol->Size = (float) Curves[i].Get_symbol_size();
				else
					myCurve->Symbol->IsVisible = false;
				myCurve->Symbol->Border->Width = (float) Curves[i].Get_line_w();
				if (Y2)
					myCurve->IsY2Axis = true;
				Curves[i].Set_npoints_plot((int) Curves[i].Get_x().size());

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
			//break;

		}

		public: void CreateGraph( ZedGraphControl ^z1 )	{
			// Get a reference to the GraphPane instance in the ZedGraphControl
			GraphPane ^myPane = z1->GraphPane;

			// lock thread
			while( 0 != System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 1) );

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
			return curve->Label->Text + " is " + pt.Y.ToString( "f3" ) + " units at X = " + pt.X.ToString( "f3" );
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
		void SaveCurves( System::Object ^sender, System::EventArgs ^e )
		{
			std::string str = "curves.u_g";
			this->chartobject_ptr->SaveCurvesToFile(str);
		}

		// Respond to a Zoom Event
		void MyZoomEvent( ZedGraphControl ^control, ZoomState ^oldState, ZoomState ^newState )
		{
			// Here we get notification everytime the user zooms
		}

		// update the chart with new data...
		private: void timer1_Tick(System::Object ^sender, System::EventArgs ^e )
		{
			LineItem  ^curve;
			ChartObject *chart = this->chartobject_ptr;
			if (chart == NULL) return;

			std::vector<CurveObject> Curves = chart->Get_CurvesCSV();
			size_t i;
			for (i = 0; i < chart->Get_CurvesPrevious().size(); i++)
			{
				Curves.push_back(chart->Get_CurvesPrevious()[i]);
			}
			for (i = 0; i < chart->Get_Curves().size(); i++)
			{
				Curves.push_back(chart->Get_Curves()[i]);
			}

			//lock for thread
			while( 0 != System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 1) );

			std::cout << std::endl << "Timer1_Tick" << std::endl;
			std::cout << "Offset:          " << chart->Get_ColumnOffset() << std::endl;

			{
				size_t i; 
				std::cout << "CSV curves:      " << chart->Get_CurvesCSV().size() << std::endl;
				for (i = 0; i < chart->Get_CurvesCSV().size(); i++)
				{
					std::cout << "\t" << i << "\t" << chart->Get_CurvesCSV()[i].Get_x().size() << "\t" << chart->Get_CurvesCSV()[i].Get_id() << std::endl;
				}
				std::cout << "Previous curves: " << chart->Get_CurvesPrevious().size() << std::endl;
				for (i = 0; i < chart->Get_CurvesPrevious().size(); i++)
				{
					std::cout << "\t" << i << "\t" << chart->Get_CurvesPrevious()[i].Get_x().size() << "\t" << chart->Get_CurvesPrevious()[i].Get_id() << std::endl;
				}
				std::cout << "Curves:          " << chart->Get_Curves().size() << std::endl;
				for (i = 0; i < chart->Get_Curves().size(); i++)
				{
					std::cout << "\t" << i << "\t" << chart->Get_Curves()[i].Get_x().size() << "\t" << chart->Get_Curves()[i].Get_id() << std::endl;
				}
			}
			std::cout << "Zedgraph curves:      " << zg1->GraphPane->CurveList->Count << std::endl;
			//this->chartobject_ptr->Get_ncurves_changed()[0] = 1;
			//std::vector<CurveObject> & Curves = chart->Get_Curves();

			if ( (Environment::TickCount - tickStart ) > this->chartobject_ptr->Get_update_time_chart()) {
				this->chartobject_ptr->Set_all_points(true);
				if (this->chartobject_ptr->Get_ncurves_changed()[0])
				{
					DefineCurves(zg1->GraphPane, zg1->GraphPane->CurveList->Count);
					this->chartobject_ptr->Set_all_points(false);
				}
				else
					// Get the graph curves...

					for (int i = 0; i < zg1->GraphPane->CurveList->Count; i++) {
						curve =  (LineItem ^) zg1->GraphPane->CurveList[i];
						// Get the PointPairList
						IPointListEdit  ^ip = (IPointListEdit^) curve->Points;
						if ((size_t) ip->Count < Curves[i].Get_x().size())
						{
							for ( size_t i2 = ip->Count; i2 < Curves[i].Get_x().size(); i2++ )
							{
								if ((LogX || LogY || LogY2) && (Curves[i].Get_x()[i2] <=0 
									|| Curves[i].Get_y()[i2] <=0))
									continue;
								else
									ip->Add(Curves[i].Get_x()[i2], Curves[i].Get_y()[i2] );
							}
						}
						//if (Curves[i].Get_npoints_plot() != Curves[i].Get_x().size())
						//	chart->Set_all_points(false);
						//else
						//	chart->Set_all_points(true);

						//for ( int i2 = Curves[i].Get_npoints_plot(); i2 < (int) Curves[i].Get_x().size(); i2++ )
						//{
						//	if ((LogX || LogY || LogY2) && (Curves[i].Get_x()[i2] <=0
						//		|| Curves[i].Get_y()[i2] <=0))
						//		continue;
						//	else
						//		ip->Add(Curves[i].Get_x()[i2], Curves[i].Get_y()[i2] );
						//}
						//Curves[i].Set_npoints_plot(Curves[i].Get_x().size());
					}
				/* explicitly reset the max in case of log scale, zedgraphs doesn't do this... */
				if ((fabs(chart->Get_axis_scale_x()[1] - NA) < 1e-3) && zg1->GraphPane->XAxis->Type == AxisType::Log)
				{
					double max = -1e99;
					for  (int i = 0; i < zg1->GraphPane->CurveList->Count; i++)
					{
						if (Curves[i].Get_x()[Curves[i].Get_x().size() - 1] > max)
							max = Curves[i].Get_x()[Curves[i].Get_x().size() - 1];
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
						if (Curves[i].Get_y()[Curves[i].Get_y().size() - 1] > max)
							max = Curves[i].Get_y()[Curves[i].Get_y().size() - 1];
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
						if (Curves[i].Get_y()[Curves[i].Get_y().size() - 1] > max)
							max = Curves[i].Get_y()[Curves[i].Get_y().size() - 1];
					}
					max += pow(10.0, log10(max / 3));
					zg1->GraphPane->Y2Axis->Scale->Max = max;
				}

				zg1->AxisChange();
				zg1->Refresh();
				tickStart = Environment::TickCount;
			}
			if (chart->Get_end_timer() && chart->Get_all_points())
			{
				timer1->Stop();
				//SaveCurvesToFile("c:\\temp\\cv.ug1");
			}

			//unlock thread
			System::Threading::Interlocked::Exchange(this->chartobject_ptr->usingResource, 0);
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
			//P_INSTANCE_POINTER DeleteCurves(); /* perhaps not even needed... */
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
			this->TopMost = true;
			System::ComponentModel::ComponentResourceManager^  resources = (gcnew System::ComponentModel::ComponentResourceManager(Form1::typeid));
			try
			{
				this->Icon = (cli::safe_cast<System::Drawing::Icon^  >(resources->GetObject(L"$this.Icon")));
			}
			catch (...)
			{
			}

			this->Load += gcnew System::EventHandler(this, &Form1::Form1_Load);
			this->Resize += gcnew System::EventHandler(this, &Form1::Form1_Resize);
			this->ResumeLayout(false);
		}
#pragma endregion
	

};
}